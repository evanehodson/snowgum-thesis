import geopandas as gpd
import numpy as np
import rasterio
import shapely.affinity as affinity
from shapely.geometry import MultiLineString


def process_australia_interior_boundaries(
    input_states_shapefile,
    reference_geotiff,
    output_svg,
    min_area_km2=150,
    pixel_buffer=100,
):
    """Processes a state/territory shapefile to extract ONLY interior state boundaries

    (excluding all external coastlines and removed island footprints), then exports
    them as a SOLID, continuous SVG perfectly aligned to the reference GeoTIFF's pixel grid.
    """
    print("\n--- Processing Interior Boundaries ---")

    # 1. Load the reference GeoTIFF to get exact pixel & spatial bounds
    print("Reading reference GeoTIFF spatial metadata...")
    with rasterio.open(reference_geotiff) as src:
        raster_pixel_width = src.width
        raster_pixel_height = src.height
        raster_min_x, raster_min_y, raster_max_x, raster_max_y = src.bounds
        raster_crs = src.crs

    # 2. Load and project the state boundaries shapefile
    print("Loading state boundaries shapefile...")
    gdf_states = gpd.read_file(input_states_shapefile)

    print(f"Projecting to match GeoTIFF CRS ({raster_crs.to_string()})...")
    gdf_states = gdf_states.to_crs(raster_crs)

    # 3. Isolate the exact same mainlands/islands used in the raster step
    print("Isolating valid mainland and large island geometries...")
    dissolved = gdf_states.dissolve()
    exploded = dissolved.explode(index_parts=False)

    # Filter landmasses using the exact same area criteria (150 sq km)
    min_area_meters = min_area_km2 * 1000 * 1000
    valid_landmasses = exploded[exploded.geometry.area >= min_area_meters]

    # Combine valid landmasses back into a single geometry mask
    landmass_mask = valid_landmasses.unary_union

    # 4. Filter original state polygons to clean up dropped islands
    print("Clipping away small island geometries from states...")
    gdf_states["geometry"] = gdf_states.geometry.intersection(landmass_mask)
    gdf_states = gdf_states[~gdf_states.geometry.is_empty]

    # 5. Extract only internal boundary lines (Solid & Unfragmented)
    print("Extracting interior-only boundary lines...")
    state_boundaries = gdf_states.geometry.boundary
    all_lines = state_boundaries.unary_union

    exterior_coastline = landmass_mask.boundary

    # FIX: Using a fraction of a meter (0.01) ensures we don't aggressively eat 
    # into the interior lines at coastal intersections, preventing fragmentation.
    interior_lines = all_lines.difference(exterior_coastline.buffer(0.01))

    if interior_lines.is_empty:
        raise ValueError(
            "No interior boundaries found! Check your shapefile attributes or geometry."
        )

    # 6. Translate and scale geometry to match pixel-space (Photoshop coordinate rules)
    print("Translating and scaling vector lines to match SVG pixel space...")

    # Step A: Translate to spatial (0,0) based on the GeoTIFF's exact bounding box
    lines_translated = affinity.translate(
        interior_lines, xoff=-raster_min_x, yoff=-raster_min_y
    )

    # Step B: Scale to canvas pixel dimensions
    scale_x = raster_pixel_width / (raster_max_x - raster_min_x)
    scale_y = raster_pixel_height / (raster_max_y - raster_min_y)
    lines_scaled = affinity.scale(
        lines_translated, xfact=scale_x, yfact=scale_y, origin=(0, 0)
    )

    # Step C: Flip Y-axis to match SVG/Photoshop's top-left origin rule
    pixel_lines = affinity.scale(
        lines_scaled, xfact=1.0, yfact=-1.0, origin=(0, 0)
    )
    pixel_lines = affinity.translate(pixel_lines, yoff=raster_pixel_height)

    # 7. Convert to SVG string data
    print("Generating perfectly aligned SVG file...")
    svg_path_data = pixel_lines.svg(scale_factor=1.0)

    # FIX: Removed stroke-dasharray="4,4" to ensure lines are completely solid
    svg_content = f"""<svg xmlns="http://www.w3.org/2000/svg" 
     width="{raster_pixel_width}" 
     height="{raster_pixel_height}" 
     viewBox="0 0 {raster_pixel_width} {raster_pixel_height}">
    <g id="state_interior_boundaries" stroke="gray" stroke-width="1.5" fill="none">
        {svg_path_data}
    </g>
</svg>"""

    with open(output_svg, "w") as f:
        f.write(svg_content)

    print(
        f"🚀 Success! Aligned SOLID interior state boundaries exported to '{output_svg}'"
    )


# --- Execution Example ---
if __name__ == "__main__":
    reference_tif = "australia_1km_filtered.tif"
    input_states_shp = (
        "STE_2021_AUST_SHP_GDA2020/STE_2021_AUST_GDA2020.shp"
    )
    output_border_svg = "australia_state_borders_aligned.svg"

    process_australia_interior_boundaries(
        input_states_shapefile=input_states_shp,
        reference_geotiff=reference_tif,
        output_svg=output_border_svg,
        min_area_km2=150,
    )