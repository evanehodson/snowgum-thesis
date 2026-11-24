"""
view_grid.py

Plots the generated grid points CSV along with the original occurrence 
and sample points in contrasting colors for visual inspection.

Usage:
    python view_grid.py \
        --input grid_points.csv \
        --occurrences ../ALA_occurences/niph_range/niph_range.csv ../ALA_occurences/pauci_range/pauci_range.csv \
        --samples eucalyptus_samples_with_elevation.csv
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt

def parse_coords(df):
    """Flexible coordinate column identifier for occurrences and samples."""
    lon_col = next((c for c in df.columns if c.lower() in ['decimallongitude', 'longitude', 'long', 'lon', 'x']), None)
    lat_col = next((c for c in df.columns if c.lower() in ['decimallatitude', 'latitude', 'lat', 'y']), None)
    
    if not lon_col or not lat_col:
        raise ValueError(f"Could not automatically identify longitude/latitude columns in dataframe. Columns found: {list(df.columns)}")
    
    return df[[lon_col, lat_col]].dropna().drop_duplicates().rename(columns={lon_col: 'lon', lat_col: 'lat'})

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", default="grid_points.csv", help="Path to the grid points CSV file")
    p.add_argument("--occurrences", nargs='*', help="Occurrence CSV files used in the grid generation")
    p.add_argument("--samples", nargs='*', help="Sample CSV files used in the grid generation")
    args = p.parse_args()

    print(f"Loading grid points from {args.input}...")
    df_grid = pd.read_csv(args.input)

    print(f"Plotting {len(df_grid):,} grid points and input datasets...")
    fig, ax = plt.subplots(figsize=(10, 10))

    # 1. Plot background grid points
    ax.scatter(
        df_grid['lon'], df_grid['lat'], 
        s=1, c='#aec7e8', alpha=0.6, 
        label=f'Grid Points (n={len(df_grid):,})'
    )

    # 2. Plot occurrence points if provided
    if args.occurrences:
        occ_dfs = []
        for file in args.occurrences:
            occ_df = pd.read_csv(file)
            occ_dfs.append(parse_coords(occ_df))
        df_occ = pd.concat(occ_dfs).dropna().drop_duplicates()
        
        ax.scatter(
            df_occ['lon'], df_occ['lat'], 
            s=8, c='#d62728', alpha=0.8, marker='o',
            label=f'Occurrences (n={len(df_occ):,})'
        )

    # 3. Plot sample points if provided
    if args.samples:
        samp_dfs = []
        for file in args.samples:
            samp_df = pd.read_csv(file)
            samp_dfs.append(parse_coords(samp_df))
        df_samp = pd.concat(samp_dfs).dropna().drop_duplicates()
        
        ax.scatter(
            df_samp['lon'], df_samp['lat'], 
            s=15, c='#2ca02c', alpha=0.9, marker='^',
            label=f'Samples (n={len(df_samp):,})'
        )
    
    ax.set_title("Grid & Input Points Inspection Map", fontsize=14, fontweight='bold')
    ax.set_xlabel("Longitude (°E)", fontsize=12)
    ax.set_ylabel("Latitude (°S)", fontsize=12)
    ax.legend(loc='upper right', frameon=True, markerscale=2)
    ax.grid(True, linestyle='--', alpha=0.4)
    
    plt.tight_layout()
    output_img = args.input.replace(".csv", "_preview.png")
    plt.savefig(output_img, dpi=300, bbox_inches='tight')
    print(f"Saved preview map to {output_img}")
    
    plt.show()

if __name__ == "__main__":
    main()