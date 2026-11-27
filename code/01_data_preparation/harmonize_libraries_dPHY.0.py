"""
Purpose:
    Standardize normalized per-tissue CSVs (`target_id` -> `gene_id`, strip trailing version
    suffixes `.1`/`.2`, standardize column to `log10MeanTPM` + placeholder `target`).

Manuscript (Peleke et al., 2026):
    Intermediate formatting step between `01_built_libraries_log10Mean_dPHY.0.py` and
    `make_labels_dPHYlogic_meanTPM.0.py` across the 13 species.

Usage & I/O:
    python code/01_data_preparation/harmonize_libraries_dPHY.0.py
    - Reads and overwrites `<BASE_RESULTS_DIR>/<species>/normalized/*_log10_{mean,median,max}_tpm.csv`

How to change or extend:
    Update `BASE_RESULTS_DIR`, `NORMALIZED_SUBDIR_NAME`, or `polish_gene_id()` regex to adapt
    to different gene identifier naming conventions.
"""
import os
import pandas as pd
import numpy as np # Still needed for NaN comparison if any, but target is "NA"
import re

# --- Configuration ---
# Base directory for input/output (where species-specific 'normalized' folders are)
BASE_RESULTS_DIR = "../data/results/"
# Name of the subdirectory containing normalized files (will be overwritten)
NORMALIZED_SUBDIR_NAME = "normalized"

# Original column name for target IDs (as read from input CSVs)
ORIGINAL_ID_COL = "target_id"
# New column name for gene IDs
NEW_GENE_ID_COL = "gene_id"
# New column name for the target classification
TARGET_CLASSIFICATION_COL = "target"
# Value to insert into the target column
TARGET_VALUE = "NA"

# Regex to find the normalization type (mean, median, max) from input column names
# Example input column: leaf_log10_median_tpm or all_tissues_log10_max_tpm
AGG_TYPE_PATTERN = re.compile(r"_log10_(mean|median|max)_tpm$", re.IGNORECASE)

def polish_gene_id(gene_id):
    """
    Removes trailing version numbers (e.g., .1, .2) from a gene ID string.
    Example: 'GENEID.1' becomes 'GENEID'.
    'GENEID.V2.1' becomes 'GENEID.V2'
    """
    if isinstance(gene_id, str):
        return re.sub(r'\.\d+$', '', gene_id)
    return gene_id

def get_new_tpm_col_name(original_value_col_name):
    """
    Determines the new TPM column name based on the original.
    e.g., 'leaf_log10_median_tpm' -> 'log10MedianTPM'
    """
    match = AGG_TYPE_PATTERN.search(original_value_col_name)
    if match:
        agg_type = match.group(1).capitalize() # Mean, Median, Max
        return f"log10{agg_type}TPM"
    else:
        # Fallback if pattern doesn't match, though it should based on previous script
        print(f"Warning: Could not parse aggregation type from '{original_value_col_name}'. Using generic name 'log10NormalizedTPM'.")
        return "log10NormalizedTPM"


def process_and_overwrite_normalized_file(filepath, species_name):
    """
    Processes a single normalized TPM file: renames columns, cleans gene IDs,
    adds a target column, and overwrites the original file.

    Args:
        filepath (str): Path to the input normalized CSV file (will be overwritten).
        species_name (str): Name of the species.

    Returns:
        bool: True if successful, False otherwise.
    """
    try:
        df = pd.read_csv(filepath)

        # Ensure the DataFrame has the expected structure (at least two columns)
        if df.shape[1] < 2:
            print(f"Warning: File {filepath} has fewer than 2 columns. Skipping.")
            return False
        
        # Identify original value column (should be the second column)
        original_value_col = df.columns[1]

        # 1. Rename target_id to gene_id
        if ORIGINAL_ID_COL not in df.columns:
            print(f"Warning: Expected column '{ORIGINAL_ID_COL}' not found in {filepath}. Assuming first column is ID. Skipping.")
            return False
        df.rename(columns={ORIGINAL_ID_COL: NEW_GENE_ID_COL}, inplace=True)

        # 2. Polish gene_id values
        df[NEW_GENE_ID_COL] = df[NEW_GENE_ID_COL].apply(polish_gene_id)

        # 3. Determine and rename the normalized TPM value column
        new_tpm_col_name = get_new_tpm_col_name(original_value_col)
        df.rename(columns={original_value_col: new_tpm_col_name}, inplace=True)

        # 4. Add the 'target' column, initialized with "NA"
        df[TARGET_CLASSIFICATION_COL] = TARGET_VALUE

        # 5. Reorder columns
        df = df[[NEW_GENE_ID_COL, new_tpm_col_name, TARGET_CLASSIFICATION_COL]]
        
        # Overwrite the original file
        try:
            df.to_csv(filepath, index=False)
            print(f"    Successfully polished and overwrote: {filepath}")
            return True
        except Exception as e:
            print(f"    Error overwriting file {filepath}: {e}")
            return False

    except FileNotFoundError:
        print(f"Error: File not found: {filepath}. Skipping.")
        return False
    except pd.errors.EmptyDataError:
        print(f"Warning: File is empty: {filepath}. Skipping.")
        return False
    except Exception as e:
        print(f"Error processing file {filepath}: {e}. Skipping.")
        return False

def main():
    """
    Main function to orchestrate the polishing of all normalized TPM files,
    overwriting them in place.
    """
    print("Starting TPM data polishing script (overwriting existing files)...")
    print(f"Looking for 'normalized' subdirectories in: {os.path.abspath(BASE_RESULTS_DIR)}")

    try:
        species_dirs = [
            d for d in os.listdir(BASE_RESULTS_DIR)
            if os.path.isdir(os.path.join(BASE_RESULTS_DIR, d))
        ]
    except FileNotFoundError:
        print(f"Error: Base results directory not found: {BASE_RESULTS_DIR}. Ensure previous scripts ran successfully.")
        return

    if not species_dirs:
        print("No species directories found in results. Exiting.")
        return

    for species_name in species_dirs:
        print(f"\nProcessing species: {species_name}")
        normalized_dir_path = os.path.join(BASE_RESULTS_DIR, species_name, NORMALIZED_SUBDIR_NAME)

        if not os.path.isdir(normalized_dir_path):
            print(f"  Normalized directory not found for {species_name}: {normalized_dir_path}. Skipping.")
            continue

        # Iterate through files in the 'normalized' directory
        for filename in os.listdir(normalized_dir_path):
            # Process files matching the expected patterns from previous normalization script
            if filename.endswith(("_log10_mean_tpm.csv", "_log10_median_tpm.csv", "_log10_max_tpm.csv")):
                filepath = os.path.join(normalized_dir_path, filename)
                print(f"  Processing file: {filename}")
                
                process_and_overwrite_normalized_file(filepath, species_name)
            else:
                # print(f"  Skipping file (doesn't match pattern): {filename}") # Optional: for debugging
                pass

    print("\nPolishing script finished.")

if __name__ == "__main__":
    main()
