"""
extractor.py
Multi-threaded extraction script for envXtractor.
Features a global authentication environment and optimized network sampling.

Usage:
    python extractor.py extraction_config.json
"""

import os
import sys
import json
import tempfile
import threading
import numpy as np
import pandas as pd
import requests
import rasterio
from rasterio.env import Env
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from dotenv import load_dotenv, find_dotenv
from tqdm import tqdm

# Load environment completely silently
load_dotenv(find_dotenv())

# Thread lock to guarantee pristine, readable console warning logs during tqdm execution
log_lock = threading.Lock()

def warn_msg(msg):
    with log_lock:
        tqdm.write(msg)

def apply_scale_offset(value, scale, offset):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan
    return float(value) * float(scale) + float(offset)


# ---------------------------------------------------------------------------
# ARCGIS FEATURE EXTRACTION METRICS
# ---------------------------------------------------------------------------
ARCGIS_METRICS = {
    "fire_count":           lambda df: len(df),
    "fire_min_interval":     lambda df: _fire_min(df),
    "fire_median_interval": lambda df: _fire_median(df),
}

def _fire_intervals(df):
    if len(df) < 2:
        return None
    df = df.sort_values("ignition_date").drop_duplicates(subset=["ignition_date"])
    return df["ignition_date"].diff().dt.days.dropna() / 365.25

def _fire_min(df):
    iv = _fire_intervals(df)
    return float(iv.min()) if iv is not None and len(iv) else np.nan

def _fire_median(df):
    iv = _fire_intervals(df)
    return float(iv.median()) if iv is not None and len(iv) else np.nan


def _query_arcgis(url, lon, lat, timeframe):
    start_year, end_year = map(int, timeframe.split("-"))
    where = f"ignition_date >= DATE '{start_year}-01-01' AND ignition_date <= DATE '{end_year}-12-31'"
    point_geometry = {"x": float(lon), "y": float(lat), "spatialReference": {"wkid": 4326}}

    params = {
        "where": where,
        "geometry": json.dumps(point_geometry),
        "geometryType": "esriGeometryPoint",
        "inSR": "4326",
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": "fire_id,ignition_date",
        "returnGeometry": "false",
        "f": "geojson",
    }
    try:
        r = requests.get(url + "/query", params=params, timeout=30)
        r.raise_for_status()
        data = r.json()
        if "error" in data:
            return pd.DataFrame()
        features = data.get("features", [])
        if not features:
            return pd.DataFrame()
        fire_df = pd.json_normalize([f["properties"] for f in features])
        fire_df = fire_df[fire_df["ignition_date"].notna()].copy()
        fire_df["ignition_date"] = pd.to_datetime(fire_df["ignition_date"], unit="ms")
        return fire_df
    except Exception:
        return None


# ---------------------------------------------------------------------------
# WORKER THREAD EXTRACTION LOGIC
# ---------------------------------------------------------------------------
def process_layer(layer, coords, total_coords, gdal_options):
    lid = layer["layer_id"]
    url = layer["url"]
    access_type = layer["access_type"]
    scale = float(layer.get("scale") or 1.0)
    offset = float(layer.get("offset") or 0.0)
    batch_size = 250  # Increased batch size for faster chunking over high-latency links

    layer_results = []

    if access_type == "cog":
        try:
            with Env(**gdal_options):
                # Major optimization: Open the file ONCE per thread rather than per batch loop
                with rasterio.open(f"/vsicurl/{url}") as src:
                    nodata = src.nodata
                    for i in range(0, total_coords, batch_size):
                        batch_coords = coords[i:i + batch_size]
                        sampled = list(src.sample(batch_coords))
                        for v in sampled:
                            raw = v[0]
                            if nodata is not None and raw == nodata:
                                layer_results.append(np.nan)
                            else:
                                layer_results.append(apply_scale_offset(raw, scale, offset))
            
            # Check if layer resulted entirely in NaN values
            if np.isnan(layer_results).all():
                warn_msg(f"⚠️ [WARNING] Layer '{lid}' processed, but returned completely NaN outputs.")
                
            return lid, layer_results

        except Exception as layer_error:
            warn_msg(f"❌ [FAIL] Layer '{lid}' network/auth error: {layer_error}")
            shortfall = total_coords - len(layer_results)
            if shortfall > 0:
                layer_results.extend([np.nan] * shortfall)
            return lid, layer_results

    elif access_type == "arcgis_feature":
        timeframe = layer.get("timeframe", "1980-2023")
        variable_id = layer.get("variable_id", "fire_count")
        
        try:
            for lon, lat in coords:
                fire_df = _query_arcgis(url, lon, lat, timeframe)
                fn = ARCGIS_METRICS.get(variable_id)
                
                if fn is None or fire_df is None:
                    layer_results.append(np.nan)
                elif fire_df.empty:
                    layer_results.append(0 if variable_id == "fire_count" else np.nan)
                else:
                    layer_results.append(fn(fire_df))
                    
            if np.isnan(layer_results).all():
                warn_msg(f"⚠️ [WARNING] ArcGIS layer '{lid}' processed, but returned entirely NaN values.")
                
            return lid, layer_results
        except Exception as arcgis_error:
            warn_msg(f"❌ [FAIL] ArcGIS layer failed: {lid} -> {arcgis_error}")
            shortfall = total_coords - len(layer_results)
            if shortfall > 0:
                layer_results.extend([np.nan] * shortfall)
            return lid, layer_results
            
    else:
        warn_msg(f"⚠️ [WARN] Skipping unknown access type for layer: {lid}")
        return lid, [np.nan] * total_coords


# ---------------------------------------------------------------------------
# MAIN POOL COORDINATOR
# ---------------------------------------------------------------------------
def extract(points_df, layers, lat_col, lon_col, max_workers=8):
    coords = list(zip(points_df[lon_col], points_df[lat_col]))
    total_coords = len(coords)
    
    extracted_columns = {}
    tern_api_key = os.environ.get("TERN_API_KEY")

    header_file_path = os.path.join(tempfile.gettempdir(), "tern_gdal_global_hdr.txt")
    with open(header_file_path, 'w') as f:
        f.write(f"X-Api-Key: {tern_api_key}\n")

    gdal_options = {
        "GDAL_HTTP_UNSAFESSL": "True",
        "GDAL_HTTP_MAX_RETRY": "5",
        "GDAL_HTTP_RETRY_DELAY": "3",
        "GDAL_HTTP_HEADER_FILE": header_file_path
    }

    try:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(process_layer, layer, coords, total_coords, gdal_options): layer["layer_id"]
                for layer in layers
            }
            
            # Use tqdm context for a clean single-line progress indicator
            with tqdm(total=len(layers), desc="Extracting Environmental Data", unit="layer") as pbar:
                for future in as_completed(futures):
                    lid, layer_results = future.result()
                    extracted_columns[lid] = layer_results
                    pbar.update(1)

    finally:
        if os.path.exists(header_file_path):
            try:
                os.remove(header_file_path)
            except Exception:
                pass

    if extracted_columns:
        extracted_df = pd.DataFrame(extracted_columns, index=points_df.index)
        overlapping_cols = [c for c in extracted_df.columns if c in points_df.columns]
        if overlapping_cols:
            points_df = points_df.drop(columns=overlapping_cols)
            
        return pd.concat([points_df, extracted_df], axis=1)
        
    return points_df


# ---------------------------------------------------------------------------
# CLI ENTRY POINT
# ---------------------------------------------------------------------------
def main():
    if not os.environ.get("TERN_API_KEY"):
        sys.exit("❌ CRITICAL ERROR: TERN_API_KEY is not loaded. Terminating immediately.")

    if len(sys.argv) < 2:
        sys.exit("Usage: python extractor.py extraction_config.json")

    config_path = Path(sys.argv[1])
    if not config_path.exists():
        sys.exit(f"Error: config file not found: {config_path}")

    with open(config_path) as f:
        config = json.load(f)

    points_file = config["points_file"]
    lat_col     = config["lat_col"]
    lon_col     = config["lon_col"]
    out_file    = config["output_file"]
    layers      = config["layers"]

    points_path = config_path.parent / points_file
    if not points_path.exists():
        sys.exit(f"Error: points file not found: {points_path}")

    points_df = pd.read_csv(points_path).dropna(subset=[lat_col, lon_col])
    points_df[lat_col] = points_df[lat_col].astype(float)
    points_df[lon_col] = points_df[lon_col].astype(float)

    result = extract(points_df, layers, lat_col, lon_col)

    out_path = config_path.parent / out_file
    result.to_csv(out_path, index=False)


if __name__ == "__main__":
    main()