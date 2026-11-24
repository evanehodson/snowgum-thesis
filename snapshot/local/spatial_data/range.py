import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from shapely.geometry import Polygon, MultiPolygon
import shapely.affinity as affinity
import rasterio

# =========================================================================
# CONFIGURATION: SPECIFY YOUR FILE PATHS HERE
# =========================================================================
CSV_PATH = "ALA_occurences/pauci_range/pauci_range.csv"
GEOTIFF_PATH = "australia_1km_filtered.tif"  # Replace with your actual GeoTIFF filename
OUTPUT_SVG_PATH = "ALA_occurences/pauci_range/pauci_range.svg"
THRESHOLD_PERCENTILE = 1  # Lower (1-2) makes it larger/rounder. Higher (5-10) hugs tighter.

# SMOOTHING FACTOR: Increase this (e.g., 5 to 15) for an even smoother, 
# more rounded "liquid" look. Decrease it to capture fine detailed curves.
SMOOTHING_BUFFER = 8 
# =========================================================================

print("1. Reading GeoTIFF spatial metadata...")
with rasterio.open(GEOTIFF_PATH) as src:
    raster_pixel_width = src.width
    raster_pixel_height = src.height
    raster_min_x, raster_min_y, raster_max_x, raster_max_y = src.bounds
    raster_crs = src.crs

print("2. Processing CSV and aligning projections...")
df = pd.read_csv(CSV_PATH).dropna(subset=["decimalLatitude", "decimalLongitude"])
geometry = gpd.points_from_xy(df["decimalLongitude"], df["decimalLatitude"])
gdf_points = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326").to_crs(raster_crs)

x = gdf_points.geometry.x.values
y = gdf_points.geometry.y.values

print("3. Computing KDE math...")
coords = np.vstack([x, y])
kde = gaussian_kde(coords)

# FIX: Dynamically scale the evaluation grid so it samples tightly.
# We sample at a much higher resolution (e.g., max 800 steps) to capture smooth geometry
grid_resolution = min(800, max(raster_pixel_width, raster_pixel_height))
X, Y = np.mgrid[raster_min_x:raster_max_x:complex(0, grid_resolution), 
                raster_min_y:raster_max_y:complex(0, grid_resolution)]

positions = np.vstack([X.ravel(), Y.ravel()])
Z = np.reshape(kde(positions).T, X.shape)

point_densities = kde(coords)
kde_threshold = np.percentile(point_densities, THRESHOLD_PERCENTILE)

# 4. Extract contour lines
fig_dummy, ax_dummy = plt.subplots()
contours = ax_dummy.contour(X, Y, Z, levels=[kde_threshold])
plt.close(fig_dummy)

polygons = []
contour_set = contours.allsegs[0] if hasattr(contours, 'allsegs') else [path.vertices for path in contours.get_paths()]
for seg in contour_set:
    if len(seg) >= 3:
        polygons.append(Polygon(seg))

if len(polygons) == 1:
    geo_shape = polygons[0]
elif len(polygons) > 1:
    geo_shape = MultiPolygon(polygons)
else:
    raise ValueError("Could not extract geometry.")

print("4. Translating and scaling to Photoshop pixel-space...")
# Step A: Translate to spatial (0,0)
shape_translated = affinity.translate(geo_shape, xoff=-raster_min_x, yoff=-raster_min_y)

# Step B: Scale to canvas pixels
scale_x = raster_pixel_width / (raster_max_x - raster_min_x)
scale_y = raster_pixel_height / (raster_max_y - raster_min_y)
shape_scaled = affinity.scale(shape_translated, xfact=scale_x, yfact=scale_y, origin=(0, 0))

# Step C: Flip Y-axis for Photoshop space
pixel_shape = affinity.scale(shape_scaled, xfact=1.0, yfact=-1.0, origin=(0, 0))
pixel_shape = affinity.translate(pixel_shape, yoff=raster_pixel_height)

# FIX: Apply a tiny mathematical vector buffer expansion and contraction.
# This forces Shapely to recalculate the sharp vertices into smooth, round Bezier splines.
pixel_shape = pixel_shape.buffer(SMOOTHING_BUFFER, join_style=1).buffer(-SMOOTHING_BUFFER, join_style=1)

print("5. Generating perfectly smooth SVG file...")
svg_path_data = pixel_shape.svg(scale_factor=1.0, fill_color="none")

svg_content = f"""<svg xmlns="http://www.w3.org/2000/svg" 
     width="{raster_pixel_width}" 
     height="{raster_pixel_height}" 
     viewBox="0 0 {raster_pixel_width} {raster_pixel_height}">
    <g id="organic_boundary_layer" stroke="purple" stroke-width="3" fill="purple" fill-opacity="0.25">
        {svg_path_data}
    </g>
</svg>"""

with open(OUTPUT_SVG_PATH, "w") as f:
    f.write(svg_content)

print(f"\n🚀 Success! Smooth, infinite-res organic SVG exported to '{OUTPUT_SVG_PATH}'")