"""
alpha_grid.py

Builds a regular lon/lat grid of points aligned to a reference raster's pixel grid.
Dynamically generates a single, unified smooth alpha shape (concave hull) from occurrence 
CSV files + mandatory sample points, and uses Natural Earth land polygons to drop water-only pixels.

Usage:
    python alpha_grid.py \
        --reference_cog https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bio12/1981-2010/CHELSA_bio12_1981-2010_V.2.1.tif \
        --occurrences sp1.csv sp2.csv \
        --samples eucalyptus_samples_with_elevation.csv \
        --output grid_points.csv
"""

import argparse
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import from_bounds
import geopandas as gpd
from shapely.geometry import Point, box, Polygon, MultiPolygon
import shapely.geometry as sg
from sklearn.cluster import DBSCAN
from scipy.spatial import Delaunay
from shapely.ops import unary_union, polygonize

def parse_coords(df):
    """Flexible coordinate column identifier for occurrences and samples."""
    lon_col = next((c for c in df.columns if c.lower() in ['decimallongitude', 'longitude', 'long', 'lon', 'x']), None)
    lat_col = next((c for c in df.columns if c.lower() in ['decimallatitude', 'latitude', 'lat', 'y']), None)
    
    if not lon_col or not lat_col:
        raise ValueError(f"Could not automatically identify longitude/latitude columns in dataframe. Columns found: {list(df.columns)}")
    
    return df[[lon_col, lat_col]].dropna().drop_duplicates().rename(columns={lon_col: 'lon', lat_col: 'lat'})

def compute_alpha_shape(points, alpha):
    """Computes a concave hull (alpha shape) from a set of xy coordinates."""
    if len(points) < 4:
        return sg.MultiPoint(points).convex_hull

    tri = Delaunay(points)
    edges = set()
    for ia, ib, ic in tri.simplices:
        pa, pb, pc = points[ia], points[ib], points[ic]
        a = np.hypot(pa[0]-pb[0], pa[1]-pb[1])
        b = np.hypot(pb[0]-pc[0], pb[1]-pc[1])
        c = np.hypot(pc[0]-pa[0], pc[1]-pa[1])
        s = (a + b + c) / 2.0
        area = np.sqrt(max(0, s * (s - a) * (s - b) * (s - c)))
        if area > 0:
            radius = (a * b * c) / (4.0 * area)
            if radius < 1.0 / alpha:
                edges.add((ia, ib))
                edges.add((ib, ic))
                edges.add((ic, ia))
    m = sg.MultiLineString([points[np.array(e)] for e in edges])
    polys = list(polygonize(m))
    if not polys:
        return sg.MultiPoint(points).convex_hull
    return unary_union(polys)

def generate_range_mask(occurrence_files, sample_files=None):
    """
    Loads occurrence and sample points, filters spatial noise/singletons using DBSCAN, 
    and generates ONE unified alpha shape from the cleaned combined dataset.
    """
    print(f"Loading occurrences from {len(occurrence_files)} file(s)...")
    dfs = []
    for file in occurrence_files:
        df = pd.read_csv(file)
        dfs.append(parse_coords(df))
        
    occ_coords = pd.concat(dfs).dropna().drop_duplicates().values

    # Load samples early so they can be cleaned and pooled together
    sample_gdf = None
    sample_coords = np.empty((0, 2))
    if sample_files:
        print(f"Loading mandatory sample file(s): {sample_files}")
        sample_dfs = []
        for sfile in sample_files:
            s_df = pd.read_csv(sfile)
            sample_dfs.append(parse_coords(s_df))
        sample_coords = pd.concat(sample_dfs).dropna().drop_duplicates().values

    # Combine all points first so DBSCAN evaluates the entire distribution uniformly
    all_raw_coords = np.vstack([occ_coords, sample_coords]) if len(sample_coords) > 0 else occ_coords
    print(f"Filtering spatial outliers and singletons across {len(all_raw_coords)} total points using DBSCAN...")

    # eps=0.5 degrees (~50km) and min_samples=3 ensures local groups (like northern points) are kept,
    # while true singletons/isolated stray points are dropped as noise (-1).
    clustering = DBSCAN(eps=0.5, min_samples=3).fit(all_raw_coords)
    core_mask = clustering.labels_ != -1
    clean_combined_coords = all_raw_coords[core_mask]

    print(f"Retained {len(clean_combined_coords)} points after noise/singleton removal ({len(all_raw_coords) - len(clean_combined_coords)} dropped).")

    # If samples were provided, separate out the surviving sample points for land-check protection
    if sample_files and len(sample_coords) > 0:
        # Re-build sample GDF using only the cleaned/retained sample coordinates
        # (or keep all input samples if they are strictly mandatory)
        sample_points = [Point(xy) for xy in sample_coords]
        sample_gdf = gpd.GeoDataFrame(geometry=sample_points, crs="EPSG:4326")

    # 3. Compute single unified alpha shape on the fully cleaned point set
    print("Generating unified alpha shape mask...")
    raw_alpha = compute_alpha_shape(clean_combined_coords, alpha=1.1)
    
    # Smooth with buffer / negative-buffer
    final_mask = raw_alpha.buffer(0.3, resolution=16).buffer(-0.15, resolution=16)
    
    return final_mask, sample_gdf

def build_grid(reference_cog, mask_poly, sample_gdf=None):
    """Builds a point grid from a COG, masked by mask_poly, and filtered against 110m land polygons."""
    minx, miny, maxx, maxy = mask_poly.bounds
    print(f"Computed bounding box from unified mask: {minx:.2f}, {miny:.2f}, {maxx:.2f}, {maxy:.2f}")

    with rasterio.Env(GDAL_HTTP_UNSAFESSL="True"):
        with rasterio.open(f"/vsicurl/{reference_cog}") as src:
            window = from_bounds(minx, miny, maxx, maxy, transform=src.transform)
            window = window.round_offsets().round_lengths()
            transform = src.window_transform(window)
            height, width = int(window.height), int(window.width)
            res_x = abs(transform.a)
            res_y = abs(transform.e)

    if height <= 0 or width <= 0:
        raise ValueError("Computed grid has zero size. Check the COG CRS and extent.")

    rows, cols = np.meshgrid(np.arange(height), np.arange(width), indexing="ij")
    xs, ys = rasterio.transform.xy(transform, rows.ravel(), cols.ravel())

    df = pd.DataFrame({
        "grid_id": [f"g{i}" for i in range(len(xs))],
        "lon": xs,
        "lat": ys,
        "row": rows.ravel(),
        "col": cols.ravel(),
    })

    print(f"Generated {len(df)} initial grid cells within bounding box. Filtering against range mask...")
    
    # 1. Filter points using the single unified spatial mask
    gdf = gpd.GeoDataFrame(df, geometry=[Point(xy) for xy in zip(df.lon, df.lat)], crs="EPSG:4326")
    df = gdf[gdf.within(mask_poly)].drop(columns="geometry").reset_index(drop=True)

    print(f"Filtered down to {len(df)} candidate points. Removing water pixels using 110m land dataset...")

    # 2. Load 110m land dataset
    land_url = "https://naturalearth.s3.amazonaws.com/110m_physical/ne_110m_land.zip"
    land_gdf = gpd.read_file(land_url)
    land_union = unary_union(land_gdf.geometry.values)

    # 3. Drop pixels that don't touch land UNLESS they contain a mandatory sample point
    half_x = res_x / 2.0
    half_y = res_y / 2.0
    
    pixel_boxes = [
        box(x - half_x, y - half_y, x + half_x, y + half_y) 
        for x, y in zip(df.lon, df.lat)
    ]
    df_pixels = gpd.GeoDataFrame(df, geometry=pixel_boxes, crs="EPSG:4326")
    
    land_mask = df_pixels.intersects(land_union)

    if sample_gdf is not None:
        sample_mask = df_pixels.intersects(unary_union(sample_gdf.geometry.values))
        final_pixel_mask = land_mask | sample_mask
    else:
        final_pixel_mask = land_mask

    df_land = df_pixels[final_pixel_mask].drop(columns="geometry").reset_index(drop=True)
    
    print(f"Final point count after removing water-only pixels: {len(df_land)}")
    return df_land, transform, height, width

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reference_cog", required=True, help="COG URL to align the grid to")
    p.add_argument("--occurrences", nargs='+', required=True, help="Occurrence CSV files to build alpha shape")
    p.add_argument("--samples", nargs='*', help="Mandatory sample CSV files to merge into alpha shape")
    p.add_argument("--output", default="grid_points.csv")
    args = p.parse_args()

    # 1. Create single unified spatial mask
    mask_poly, sample_gdf = generate_range_mask(
        occurrence_files=args.occurrences, 
        sample_files=args.samples
    )

    # 2. Build and mask the grid
    df, transform, height, width = build_grid(args.reference_cog, mask_poly, sample_gdf)
    
    # 3. Save outputs
    df.to_csv(args.output, index=False)
    print(f"Wrote {len(df)} grid points to {args.output}")
    print(f"Full grid shape: {height} rows x {width} cols")

    meta_path = args.output.replace(".csv", "_meta.txt")
    with open(meta_path, "w") as f:
        f.write(f"{height},{width}\n")
        f.write(",".join(map(str, transform.to_gdal())) + "\n")
    print(f"Saved grid shape/transform to {meta_path}")

if __name__ == "__main__":
    main()