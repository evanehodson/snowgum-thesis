#!/usr/bin/env python3
"""
Download target CMIP6 annual AI/PET zip files from ScienceDB to GADI.

Target set: GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0, UKESM1-0-LL
SSPs: 126, 370, 585 (GFDL-ESM4 x SSP585 skipped, not available in source)
Period: 2041-2060, annual AI/PET product

Plus a historical baseline file (WorldClim_2_1_1970-2000).

Run on a GADI login node, or via a copyq job (copyq allows sequential
transfers only, hence no threading here).
"""

import os
import requests

# --- Config -------------------------------------------------------------
OUTPUT_DIR = "/scratch/xe2/eh4663/spatial/Global_Aridity_PET/"

# fileId + fileName pulled from the ScienceDB manifest, filtered to the
# exact target combination (model, ssp) for the annual AI/PET 2041-2060 set.
FILES = [
    ("fa33bf88d32f18d1ea8a7031cff907b1", "UKESM1-0-LL_585_2041-2060.zip"),
    ("1a6ad4e76b28d45b3a67a6bf3c7a9954", "IPSL-CM6A-LR_126_2041-2060.zip"),
    ("d3bbd27415ba7137e788f54249eb8c69", "IPSL-CM6A-LR_370_2041-2060.zip"),
    ("c8fc50fe0024a17dbcbde82c85039aa3", "IPSL-CM6A-LR_585_2041-2060.zip"),
    ("8b7587a56a3d284a13ea4d0d20ef752f", "MRI-ESM2-0_370_2041-2060.zip"),
    ("f40d739ea25faf341528e70f6457bd91", "MRI-ESM2-0_126_2041-2060.zip"),
    ("d7302c36b915b499a686efb937100ed2", "UKESM1-0-LL_370_2041-2060.zip"),
    ("4d55b8e16a3079ee2c28272a7d3ef7e8", "GFDL-ESM4_126_2041-2060.zip"),
    ("8ce4efe0e101e3a18fe1bf06a4d51ebc", "MRI-ESM2-0_585_2041-2060.zip"),
    ("9cfaaa2d5c164677cb14b6cbee19c4b3", "MPI-ESM1-2-HR_370_2041-2060.zip"),
    ("189533443bf5ffe785405ef47975344d", "UKESM1-0-LL_126_2041-2060.zip"),
    ("8f0467f0d0af98d124ee36d658613bef", "MPI-ESM1-2-HR_585_2041-2060.zip"),
    ("2fd0b6334765c37c1205d7d645149950", "GFDL-ESM4_370_2041-2060.zip"),
    ("2a2c90a85cb86f7f08036afdfad9b20f", "MPI-ESM1-2-HR_126_2041-2060.zip"),
]

BASE_PATH = "/V7/CMIP6_Models/CIMP6_PET_AI_annual_2041-2060"

# Historical baseline, different folder in the ScienceDB tree.
EXTRA_FILES = [
    (
        "dd77d8e4b0996c9620333e71b616acff",
        "WorldClim_2_1_1970-2000.zip",
        "/V7/Historical_Datasets/Historical_PET_AI_annual_1960-2000",
    ),
]


def build_url(file_id, filename, base_path=BASE_PATH):
    path = f"{base_path}/{filename}"
    return f"https://download.scidb.cn/download?fileId={file_id}&path={path}&fileName={filename}"


def download(url, filepath, chunk_size=1 << 20):
    """Stream download with basic retry-friendly behaviour."""
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        tmp_path = filepath + ".part"
        total = int(r.headers.get("content-length", 0))
        written = 0
        with open(tmp_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    written += len(chunk)
        os.rename(tmp_path, filepath)
        if total and written != total:
            print(f"   ! warning: expected {total} bytes, got {written} bytes")


def fetch_one(file_id, filename, base_path=BASE_PATH):
    filepath = os.path.join(OUTPUT_DIR, filename)

    if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
        print(f"-> Already downloaded: {filename}")
        return True

    url = build_url(file_id, filename, base_path)
    print(f"-> Downloading: {filename}")
    try:
        download(url, filepath)
        size_mb = os.path.getsize(filepath) / 1e6
        print(f"   done: {filename} ({size_mb:.1f} MB)")
        return True
    except Exception as e:
        print(f"   FAILED: {filename}: {e}")
        return False


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    all_jobs = [(fid, fname, BASE_PATH) for fid, fname in FILES] + [
        (fid, fname, base) for fid, fname, base in EXTRA_FILES
    ]
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Downloading {len(all_jobs)} target files sequentially...\n")

    failures = []
    for file_id, filename, base_path in all_jobs:
        ok = fetch_one(file_id, filename, base_path)
        if not ok:
            failures.append(filename)

    print("\nDownload run complete.")
    if failures:
        print(f"{len(failures)} file(s) failed:")
        for f in failures:
            print(f"  - {f}")
    else:
        print("All files downloaded successfully.")


if __name__ == "__main__":
    main()