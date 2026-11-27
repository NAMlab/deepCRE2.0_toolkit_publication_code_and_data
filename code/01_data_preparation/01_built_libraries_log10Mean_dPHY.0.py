"""
Purpose:
    Compute log10(mean TPM + 1) expression values per tissue (leaf, root) and across all
    tissues from each species' TPM matrix CSV.

Manuscript (Peleke et al., 2026):
    Produces the normalized expression metric plotted in Figure 1b and consumed downstream
    by `harmonize_libraries_dPHY.0.py` and `make_labels_dPHYlogic_meanTPM.0.py`.

Usage & I/O:
    python code/01_data_preparation/01_built_libraries_log10Mean_dPHY.0.py
    - Reads `<BASE_INPUT_DIR>/<species>/<species>_tpm_matrix.csv`
    - Writes `<BASE_INPUT_DIR>/<species>/normalized/<species>_<tissue>_log10_mean_tpm.csv`

How to change or extend:
    Edit `BASE_INPUT_DIR` and `NORMALIZED_OUTPUT_SUBDIR` at the top of the file, or adjust
    `normalize_log10_mean_plus_one()` to test alternative expression transformations.
"""
import os
import pandas as pd
import numpy as np
import re # For parsing column names

# --- Configuration ---
# Base directory for input data (where species-specific TPM matrices are)
# This is the output directory of the first aggregation script.
BASE_INPUT_DIR = "../data/results/"
# Suffix for the input TPM matrix files from the aggregation script
TPM_MATRIX_SUFFIX = "_tpm_matrix.csv"
# Subdirectory within each species' result folder to save normalized data
NORMALIZED_OUTPUT_SUBDIR = "normalized"
# Column name for target IDs in the input TPM matrix
TARGET_ID_COL = "target_id" # This is the column name from the aggregated matrix

def load_tpm_matrix(species_name):
    """
    Loads the TPM matrix for a given species.

    Args:
        species_name (str): The name of the species.

    Returns:
        pandas.DataFrame: The TPM matrix with target_id as index, or None if not found.
    """
    tpm_matrix_filename = f"{species_name}{TPM_MATRIX_SUFFIX}"
    # Path to the aggregated TPM matrix (e.g., ../data/results/hordeum_vulgare/hordeum_vulgare_tpm_matrix.csv)
    tpm_matrix_path = os.path.join(BASE_INPUT_DIR, species_name, tpm_matrix_filename)

    if not os.path.exists(tpm_matrix_path):
        print(f"Warning: TPM matrix not found for {species_name} at {tpm_matrix_path}. Skipping.")
        return None
    try:
        df = pd.read_csv(tpm_matrix_path)
        if TARGET_ID_COL not in df.columns:
            print(f"Warning: '{TARGET_ID_COL}' not found in {tpm_matrix_path}. Please ensure the aggregation script ran correctly. Skipping.")
            return None
        df.set_index(TARGET_ID_COL, inplace=True)
        return df
    except pd.errors.EmptyDataError:
        print(f"Warning: TPM matrix file is empty: {tpm_matrix_path}. Skipping.")
        return None
    except Exception as e:
        print(f"Error loading TPM matrix {tpm_matrix_path}: {e}. Skipping.")
        return None

def get_tissue_types_from_columns(df_columns):
    """
    Extracts unique tissue types from the DataFrame column names.
    Assumes column format like 'SRRxxxx_tpm_tissue'.
    """
    tissue_types = set()
    # Regex to capture the tissue part: SRR..._tpm_(tissue_name)
    pattern = re.compile(r"^.+_tpm_([a-zA-Z0-9_]+)$")
    for col in df_columns:
        match = pattern.match(col)
        if match:
            tissue_types.add(match.group(1))
    return list(tissue_types)

def normalize_log10_mean_plus_one(data_series):
    """
    Applies log10(value + 1) normalization to a pandas Series.
    The input data_series should already be the result of a .mean() operation.
    """
    return np.log10(data_series + 1)

def save_normalized_data(df_to_save, species_name, output_filename_stem, normalized_value_col_name):
    """
    Saves the normalized DataFrame to a CSV file.
    The DataFrame should have target_id (or the new gene_id) as index and one column for normalized values.

    Args:
        df_to_save (pandas.DataFrame): DataFrame to save (index is target_id).
        species_name (str): Name of the species.
        output_filename_stem (str): The part of the filename identifying the normalization type (e.g., tissue_name or 'all_tissues').
        normalized_value_col_name (str): The name for the column containing normalized values.
    """
    output_dir = os.path.join(BASE_INPUT_DIR, species_name, NORMALIZED_OUTPUT_SUBDIR)
    os.makedirs(output_dir, exist_ok=True)

    # Update filename to reflect log10 mean normalization
    output_filename = f"{species_name}_{output_filename_stem}_log10_mean_tpm.csv"
    output_filepath = os.path.join(output_dir, output_filename)

    # Reset index to make target_id a column for the CSV
    df_to_save_csv = df_to_save.reset_index()
    # Ensure the columns are named correctly: the ID column (still target_id at this stage) and the value column
    df_to_save_csv.columns = [TARGET_ID_COL, normalized_value_col_name]


    try:
        df_to_save_csv.to_csv(output_filepath, index=False)
        print(f"  Saved: {output_filepath}")
    except Exception as e:
        print(f"  Error saving normalized data to {output_filepath}: {e}")

def main():
    """
    Main function to orchestrate loading, normalizing (log10 mean + 1), and saving TPM data.
    """
    print("Starting TPM normalization script (log10 mean)...")
    print(f"Looking for species TPM matrices in: {os.path.abspath(BASE_INPUT_DIR)}")

    try:
        # Species are identified by subdirectories in BASE_INPUT_DIR
        # (excluding common non-species dirs like 'normalized' or 'polished_normalized' if scripts are re-run)
        species_list = [
            d for d in os.listdir(BASE_INPUT_DIR)
            if os.path.isdir(os.path.join(BASE_INPUT_DIR, d)) and d not in [NORMALIZED_OUTPUT_SUBDIR, "polished_normalized"]
        ]
    except FileNotFoundError:
        print(f"Error: Base input directory not found: {BASE_INPUT_DIR}. Ensure the TPM aggregation script ran successfully.")
        return

    if not species_list:
        print(f"No species directories found in {BASE_INPUT_DIR} (excluding '{NORMALIZED_OUTPUT_SUBDIR}' and 'polished_normalized'). Exiting.")
        return

    for species_name in species_list:
        print(f"\nProcessing species: {species_name}")
        tpm_df = load_tpm_matrix(species_name) # This loads the matrix with 'target_id' as index

        if tpm_df is None or tpm_df.empty:
            continue

        # --- 1. Per-tissue normalization (log10 mean + 1) ---
        tissue_types = get_tissue_types_from_columns(tpm_df.columns)
        if not tissue_types:
            print(f"  No tissue-specific TPM columns found for {species_name} (e.g., SRRxxxx_tpm_tissue). Skipping per-tissue normalization.")
        else:
            print(f"  Found tissues: {', '.join(tissue_types)}")

        for tissue in tissue_types:
            tissue_cols_pattern = f"_tpm_{tissue}$" # Ends with _tpm_tissue
            tissue_tpm_cols = [col for col in tpm_df.columns if re.search(tissue_cols_pattern, col)]

            if not tissue_tpm_cols:
                print(f"  Warning: No TPM columns found for tissue '{tissue}' in {species_name}. Skipping.")
                continue

            # Calculate MEAN TPM for the current tissue per target_id
            mean_tissue_tpm = tpm_df[tissue_tpm_cols].mean(axis=1)

            # Apply log10(mean + 1)
            normalized_tissue_tpm = normalize_log10_mean_plus_one(mean_tissue_tpm)
            
            df_to_save = normalized_tissue_tpm.to_frame() # Convert Series to DataFrame
            # Update output column name
            output_col_name = f"{tissue}_log10_mean_tpm"
            save_normalized_data(df_to_save, species_name, tissue, output_col_name)

        # --- 2. Overall mean normalization (across all tissues and replicates, log10 mean + 1) ---
        if tpm_df.empty or tpm_df.shape[1] == 0:
            print(f"  No TPM data columns found in the matrix for {species_name}. Skipping overall normalization.")
        else:
            # Calculate MEAN TPM across all samples for each target_id
            overall_mean_tpm = tpm_df.mean(axis=1)

            # Apply log10(mean + 1)
            normalized_overall_tpm = normalize_log10_mean_plus_one(overall_mean_tpm)

            df_to_save_overall = normalized_overall_tpm.to_frame() # Convert Series to DataFrame
            # Update output column name
            output_col_name_overall = "all_tissues_log10_mean_tpm"
            save_normalized_data(df_to_save_overall, species_name, "all_tissues", output_col_name_overall)

    print("\nNormalization script (log10 mean) finished.")

if __name__ == "__main__":
    main()
