import os
import glob
import pandas as pd

# 1. Define the root directory where your folders live
# Since merge.py is right next to the folders, we can use the current working directory.
root_dir = os.path.dirname(os.path.abspath(__file__))

# 2. Set up a mapping to convert the folder/file context into your 7 clean target taxa labels
taxon_mapping = {
    "pauci_range": "E. pauciflora subsp. pauciflora",
    "niph_range" : "niphophila",
    "deb_range"  : "debeuzivili",
    "acer_range" : "acerina",
    "hed_range"  : "hedraia",
    "pavi_range" : "parvifructa",
    "lac_range"  : "E. lacrimens"
}

# This list will temporarily hold the individual dataframes
all_dfs = []

print("Starting to compile ALA occurrence files...")

# 3. Loop through your designated folders
for folder_name, standardized_name in taxon_mapping.items():
    folder_path = os.path.join(root_dir, folder_name)
    
    if not os.path.isdir(folder_path):
        print(f"Warning: Expected folder '{folder_name}' not found. Skipping...")
        continue
        
    # Dynamically find the main occurrence file matching the pattern *_range.csv
    # This safely ignores citation.csv and headings.csv
    search_pattern = os.path.join(folder_path, f"{folder_name.split('_')[0]}_range.csv")
    matched_files = glob.glob(search_pattern)
    
    if not matched_files:
        print(f"Warning: Could not find occurrence file in {folder_name}/. Skipping...")
        continue
        
    file_path = matched_files[0]
    print(f"--> Processing: {folder_name} ({standardized_name})")
    
    try:
        # Read the csv file. 
        # Using low_memory=False because ALA files have many columns with mixed data types.
        temp_df = pd.read_csv(file_path, low_memory=False)
        
        # Inject the standardized taxon name as a dedicated tracking column
        temp_df['taxon_group'] = standardized_name
        
        # Keep only essential spatial columns + your taxon tracker to prevent memory issues 
        # (ALA sheets have ~150+ columns, most of which are empty or metadata).
        # You can add or modify columns in this list depending on what you need.
        columns_to_keep = [
            'taxon_group', 'scientificName', 'decimalLatitude', 'decimalLongitude',
            'stateProvince', 'year', 'month', 'day', 'basisOfRecord'
        ]
        
        # Intersection ensures we don't break if one specific file lacks an optional column
        available_cols = [col for col in columns_to_keep if col in temp_df.columns]
        temp_df_filtered = temp_df[available_cols]
        
        all_dfs.append(temp_df_filtered)
        
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")

# 4. Concatenate all rows into a single dataframe
if all_dfs:
    combined_df = pd.concat(all_dfs, ignore_index=True)
    
    # Optional: Clean out any records that are missing spatial coordinates
    initial_count = len(combined_df)
    combined_df = combined_df.dropna(subset=['decimalLatitude', 'decimalLongitude'])
    dropped_count = initial_count - len(combined_df)
    
    # Save the output file to your desktop folder
    output_filename = "combined_ala_occurrences.csv"
    output_path = os.path.join(root_dir, output_filename)
    combined_df.to_csv(output_path, index=False)
    
    print("\n--- Success! ---")
    print(f"Total occurrences combined: {len(combined_df)}")
    if dropped_count > 0:
        print(f"Dropped {dropped_count} rows due to missing latitude/longitude coordinates.")
    print(f"Master file saved as: {output_path}")
else:
    print("\nNo files were successfully processed.")