import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.transform import from_origin


def process_australia_shapefile(input_shapefile, output_raster, pixel_buffer=100):
    # 1. Load the shapefile
    print("Loading shapefile...")
    gdf = gpd.read_file(input_shapefile)

    # 2. Project to Australian Albers (EPSG:3577)
    print("Projecting to Australian Albers (EPSG:3577)...")
    gdf_albers = gdf.to_crs(epsg=3577)

    # 3. Clean up small islands directly in Vector space first
    print("Filtering out islands smaller than 150 sq km...")
    # 150 sq km = 150,000,000 square meters
    min_area_meters = 150 * 1000 * 1000

    # Filter out features smaller than our threshold
    # (Explode handles MultiPolygons properly so individual small islands drop out)
    gdf_exploded = gdf_albers.explode(index_parts=False)
    gdf_filtered = gdf_exploded[gdf_exploded.geometry.area >= min_area_meters]

    if gdf_filtered.empty:
        raise ValueError(
            "No geometries left after filtering! Check your area threshold."
        )

    # 4. Define raster grid properties based on the CLEANED vector bounds
    pixel_size = 1000  # meters (1km)
    xmin, ymin, xmax, ymax = gdf_filtered.total_bounds

    # 5. Apply the pixel buffer directly to the geographic bounds
    # Since resolution is 1000m, buffer meters = pixel_buffer * 1000
    buffer_meters = pixel_buffer * pixel_size
    xmin -= buffer_meters
    ymin -= buffer_meters
    xmax += buffer_meters
    ymax += buffer_meters

    # Calculate final width and height in pixels
    width = int(np.ceil((xmax - xmin) / pixel_size))
    height = int(np.ceil((ymax - ymin) / pixel_size))

    # Re-align grid perfectly to pixel dimensions
    xmax = xmin + (width * pixel_size)
    ymin = ymax - (height * pixel_size)

    # Create the correct final affine transform matrix
    transform = from_origin(xmin, ymax, pixel_size, pixel_size)

    print(
        f"Raster dimensions: {width}x{height} (Includes an unconstrained {pixel_buffer}px / {pixel_buffer}km buffer)"
    )

    # 6. Rasterize the remaining large geometries
    print("Rasterizing filtered shapefile...")
    shapes = ((geom, 1) for geom in gdf_filtered.geometry if geom is not None)
    final_array = rasterize(
        shapes=shapes,
        out_shape=(height, width),
        transform=transform,
        fill=0,
        all_touched=False,
        dtype=np.uint8,
    )

    # 7. Save the final raster to a GeoTIFF file
    print(f"Writing output to {output_raster}...")
    with rasterio.open(
        output_raster,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype=np.uint8,
        crs=gdf_filtered.crs,
        transform=transform,
        compress="lzw",
    ) as dst:
        dst.write(final_array, 1)

    print("Processing complete!")


# --- Execution Example ---
if __name__ == "__main__":
    input_shp = "../../AUS_2021_AUST_SHP_GDA2020/AUS_2021_AUST_GDA2020.shp"
    output_tif = "australia_1km_filtered.tif"

    # This will give you exactly 100 pixels (100km) of space on all 4 sides
    process_australia_shapefile(input_shp, output_tif, pixel_buffer=100)