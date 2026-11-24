#!/usr/bin/env python3
"""
Appends Aridity Index (AI) as a new predictor column to:
  - current_covariates.csv       (site-level, columns: lat, long, ...)
  - current_covariate_grid.csv   (grid-level, columns: lat, lon, ...)
  - future_covariate_grid.csv    (grid-level, columns: lat, lon, ...)

Naming convention (matches existing bio_xx columns so gdm_final's
shared_predictors intersection and process_grid_ssp regex both work):
  baseline : AI_1981-2010
  future   : AI_{period}_{model}_ssp{ssp}   e.g. AI_2041-2060_IPSL-CM6A-LR_ssp585

This OVERWRITES each input CSV in place with the new column appended.
Back up your CSVs before running if you want to keep the originals untouched.
"""

import os
import subprocess
import numpy as np
import pandas as pd
from osgeo import gdal

gdal.UseExceptions()
print = lambda *a, **k: __builtins__.print(*a, **{**k, "flush": True})

ZIP_DIR = "/scratch/xe2/eh4663/spatial/Global_Aridity_PET/"
EXTRACT_DIR = "/scratch/xe2/eh4663/spatial/Global_Aridity_PET/extracted/"

BASELINE_ZIP = "WorldClim_2_1_1970-2000.zip"
BASELINE_PERIOD_LABEL = "1981-2010"   # matches existing bio_xx_1981-2010 columns

FUTURE_PERIOD_LABEL = "2041-2060"     # dash-separated, matches bio_xx_{period}_{model}_ssp{ssp}

CMIP6_FILES = [
    ("GFDL-ESM4", "126", "GFDL-ESM4_126_2041-2060.zip"),
    ("GFDL-ESM4", "370", "GFDL-ESM4_370_2041-2060.zip"),
    ("IPSL-CM6A-LR", "126", "IPSL-CM6A-LR_126_2041-2060.zip"),
    ("IPSL-CM6A-LR", "370", "IPSL-CM6A-LR_370_2041-2060.zip"),
    ("IPSL-CM6A-LR", "585", "IPSL-CM6A-LR_585_2041-2060.zip"),
    ("MPI-ESM1-2-HR", "126", "MPI-ESM1-2-HR_126_2041-2060.zip"),
    ("MPI-ESM1-2-HR", "370", "MPI-ESM1-2-HR_370_2041-2060.zip"),
    ("MPI-ESM1-2-HR", "585", "MPI-ESM1-2-HR_585_2041-2060.zip"),
    ("MRI-ESM2-0", "126", "MRI-ESM2-0_126_2041-2060.zip"),
    ("MRI-ESM2-0", "370", "MRI-ESM2-0_370_2041-2060.zip"),
    ("MRI-ESM2-0", "585", "MRI-ESM2-0_585_2041-2060.zip"),
    ("UKESM1-0-LL", "126", "UKESM1-0-LL_126_2041-2060.zip"),
    ("UKESM1-0-LL", "370", "UKESM1-0-LL_370_2041-2060.zip"),
    ("UKESM1-0-LL", "585", "UKESM1-0-LL_585_2041-2060.zip"),
]
# Note: GFDL-ESM4 has no ssp585 entry here, carried over from your original
# CMIP6_FILES list (14, not 15 combos) - presumably matches what you have
# for the existing bio_xx future columns already. If that's not the case,
# process_grid_ssp will produce an NA-heavy AI column for that model/ssp.

TARGETS = [
    dict(kind="current",      path="current_covariates.csv",     lat_col="lat", lon_col="long"),
    dict(kind="current_grid", path="current_covariate_grid.csv", lat_col="lat", lon_col="lon"),
    dict(kind="future_grid",  path="future_covariate_grid.csv",  lat_col="lat", lon_col="lon"),
]


def ensure_extracted(zip_name):
    stem = zip_name.replace(".zip", "")
    out_dir = os.path.join(EXTRACT_DIR, stem)
    ai_path = os.path.join(out_dir, "aridity_index.tif")

    if os.path.exists(ai_path):
        return ai_path

    os.makedirs(out_dir, exist_ok=True)
    zip_path = os.path.join(ZIP_DIR, zip_name)
    print(f"  extracting {zip_name} ...")
    result = subprocess.run(
        ["unzip", "-o", zip_path, "aridity_index.tif", "-d", out_dir],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr)

    return ai_path


def sample_raster(tif_path, lats, lons):
    print(f"    opening {os.path.basename(tif_path)} ...")
    ds = gdal.Open(tif_path)
    band = ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()

    gt = ds.GetGeoTransform()
    origin_x, pixel_w, _, origin_y, _, pixel_h = gt
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize

    col = np.floor((lons - origin_x) / pixel_w).astype(int)
    row = np.floor((lats - origin_y) / pixel_h).astype(int)

    valid = (col >= 0) & (col < xsize) & (row >= 0) & (row < ysize)
    values = np.full(len(lats), np.nan, dtype=np.float64)

    if valid.any():
        col_min, col_max = col[valid].min(), col[valid].max()
        row_min, row_max = row[valid].min(), row[valid].max()
        win_xsize = col_max - col_min + 1
        win_ysize = row_max - row_min + 1

        print(f"    reading window ({win_xsize}x{win_ysize}) ...")
        window = band.ReadAsArray(int(col_min), int(row_min), int(win_xsize), int(win_ysize))

        rel_col = col[valid] - col_min
        rel_row = row[valid] - row_min
        values[valid] = window[rel_row, rel_col].astype(np.float64)

    if nodata is not None:
        values[values == nodata] = np.nan

    ds = None
    return values


def to_int(arr):
    return pd.Series(arr).round().astype("Int64")


def append_baseline_ai(target):
    df = pd.read_csv(target["path"])
    lats = df[target["lat_col"]].to_numpy()
    lons = df[target["lon_col"]].to_numpy()

    colname = f"AI_{BASELINE_PERIOD_LABEL}"
    print(f"[{target['path']}] extracting baseline AI -> {colname}")
    ai_path = ensure_extracted(BASELINE_ZIP)
    ai_vals = sample_raster(ai_path, lats, lons)

    df[colname] = to_int(ai_vals)
    df.to_csv(target["path"], index=False)
    print(f"  wrote {target['path']} with column {colname}")


def append_future_ai(target):
    df = pd.read_csv(target["path"])
    lats = df[target["lat_col"]].to_numpy()
    lons = df[target["lon_col"]].to_numpy()

    for model, ssp, zip_name in CMIP6_FILES:
        colname = f"AI_{FUTURE_PERIOD_LABEL}_{model}_ssp{ssp}"
        print(f"[{target['path']}] extracting {model} SSP{ssp} -> {colname}")
        ai_path = ensure_extracted(zip_name)
        ai_vals = sample_raster(ai_path, lats, lons)
        df[colname] = to_int(ai_vals)

    df.to_csv(target["path"], index=False)
    print(f"  wrote {target['path']} with {len(CMIP6_FILES)} AI columns")


def main():
    for target in TARGETS:
        if not os.path.exists(target["path"]):
            print(f"  skipping {target['path']} (not found in cwd)")
            continue
        if target["kind"] in ("current", "current_grid"):
            append_baseline_ai(target)
        elif target["kind"] == "future_grid":
            append_future_ai(target)


if __name__ == "__main__":
    main()