"""
Purpose:
    Aggregate per-run kallisto `abundance.tsv` files across RNA-seq experiments into one
    gene x sample TPM matrix CSV per species (`<species>_tpm_matrix.csv`).

Manuscript (Peleke et al., 2026):
    Original upstream TPM-matrix builder producing `data/<species>/<species>_tpm_matrix.csv.gz`
    (preserved with its original relative `../data/` cluster paths).

Usage & I/O:
    python code/01_data_preparation/00_built_read_count_database_dPHY.0.py
    - Reads `<BASE_DATA_DIR>/<species>/<MAPPING_OUTPUT_DIR_NAME>/<tissue>/<run>/abundance.tsv`
    - Writes `<RESULTS_DIR>/<species>/<species>_tpm_matrix.csv`

How to change or extend:
    Update `BASE_DATA_DIR`, `RESULTS_DIR`, `MAPPING_OUTPUT_DIR_NAME`, or column names
    (`TARGET_ID_COL`, `TPM_COL`) at the top of this file before running locally.
"""
import os
import pandas as pd

# --- Configuration ---
# Base directory for input data, relative to the script's location in 'bin'
BASE_DATA_DIR = "../data/"
# Base directory for output results, relative to the script's location
RESULTS_DIR = os.path.join(BASE_DATA_DIR, "results")
# Name of the mapping output directory
MAPPING_OUTPUT_DIR_NAME = "mappingoutput"
# Name of the abundance file
ABUNDANCE_FILE_NAME = "abundance.tsv"
# Columns to extract from the abundance file
TARGET_ID_COL = "target_id"
TPM_COL = "tpm"

def process_single_abundance_file(file_path, srr_accession, tissue_type):
    """
    Reads a single abundance.tsv file, extracts target_id and tpm,
    and renames the tpm column.

    Args:
        file_path (str): Path to the abundance.tsv file.
        srr_accession (str): The SRR accession ID (directory name).
        tissue_type (str): The tissue type (directory name).

    Returns:
        pandas.DataFrame: DataFrame with target_id as index and one tpm column,
                          or None if an error occurs.
    """
    try:
        # Read the abundance file
        df = pd.read_csv(file_path, sep='\t', usecols=[TARGET_ID_COL, TPM_COL])

        # Check if essential columns are present
        if TARGET_ID_COL not in df.columns or TPM_COL not in df.columns:
            print(f"Warning: Missing '{TARGET_ID_COL}' or '{TPM_COL}' in {file_path}. Skipping.")
            return None

        # Construct the new TPM column name
        new_tpm_col_name = f"{srr_accession}_tpm_{tissue_type}"
        df.rename(columns={TPM_COL: new_tpm_col_name}, inplace=True)

        # Set target_id as index
        df.set_index(TARGET_ID_COL, inplace=True)

        return df

    except FileNotFoundError:
        print(f"Warning: File not found: {file_path}. Skipping.")
        return None
    except pd.errors.EmptyDataError:
        print(f"Warning: File is empty: {file_path}. Skipping.")
        return None
    except Exception as e:
        print(f"Error processing file {file_path}: {e}. Skipping.")
        return None

def process_species(species_name):
    """
    Processes all abundance files for a given species and merges them.

    Args:
        species_name (str): Name of the species (directory name).

    Returns:
        pandas.DataFrame: A merged DataFrame containing all TPM values for the species,
                          or None if no data could be processed.
    """
    print(f"\nProcessing species: {species_name}...")
    species_path = os.path.join(BASE_DATA_DIR, species_name)
    all_sample_dfs = [] # To store DataFrames from each sample

    # Iterate through tissue types (leaf, root, shoot, etc.)
    try:
        tissue_types = [d for d in os.listdir(species_path) if os.path.isdir(os.path.join(species_path, d))]
    except FileNotFoundError:
        print(f"Warning: Species directory not found: {species_path}. Skipping species.")
        return None

    for tissue_type in tissue_types:
        tissue_path = os.path.join(species_path, tissue_type)
        mapping_output_path = os.path.join(tissue_path, MAPPING_OUTPUT_DIR_NAME)

        if not os.path.isdir(mapping_output_path):
            # print(f"Info: '{MAPPING_OUTPUT_DIR_NAME}' not found in {tissue_path}. Skipping this tissue path.")
            continue

        # Iterate through SRR accession directories
        try:
            srr_accessions = [d for d in os.listdir(mapping_output_path) if os.path.isdir(os.path.join(mapping_output_path, d))]
        except FileNotFoundError:
            print(f"Warning: Mapping output directory not found: {mapping_output_path} (expected). Skipping.")
            continue
            
        for srr_accession in srr_accessions:
            abundance_file_path = os.path.join(mapping_output_path, srr_accession, ABUNDANCE_FILE_NAME)
            
            print(f"  Reading: {abundance_file_path}")
            sample_df = process_single_abundance_file(abundance_file_path, srr_accession, tissue_type)
            
            if sample_df is not None:
                all_sample_dfs.append(sample_df)

    if not all_sample_dfs:
        print(f"No valid abundance data found for species: {species_name}.")
        return None

    # Merge all DataFrames for the current species
    # Concatenating is generally more efficient than repeated merging for many DataFrames
    print(f"Merging {len(all_sample_dfs)} samples for {species_name}...")
    try:
        merged_species_df = pd.concat(all_sample_dfs, axis=1, join='outer')
        # `axis=1` concatenates columns (samples)
        # `join='outer'` keeps all target_ids, filling missing with NaN
    except Exception as e:
        print(f"Error merging DataFrames for species {species_name}: {e}")
        return None
        
    print(f"Successfully merged data for {species_name}.")
    return merged_species_df

def main():
    """
    Main function to orchestrate the processing of all species.
    """
    print("Starting TPM data aggregation script...")

    # Create the main results directory if it doesn't exist
    os.makedirs(RESULTS_DIR, exist_ok=True)
    print(f"Output will be saved in: {os.path.abspath(RESULTS_DIR)}")

    # Get list of species directories
    try:
        all_items_in_data_dir = os.listdir(BASE_DATA_DIR)
        species_dirs = [
            item for item in all_items_in_data_dir
            if os.path.isdir(os.path.join(BASE_DATA_DIR, item)) and item != "results"
        ]
    except FileNotFoundError:
        print(f"Error: Base data directory not found: {BASE_DATA_DIR}. Please ensure the script is run from the 'bin' directory and the data structure is correct.")
        return
        
    if not species_dirs:
        print(f"No species directories found in {BASE_DATA_DIR}.")
        return

    print(f"Found species: {', '.join(species_dirs)}")

    for species_name in species_dirs:
        merged_df = process_species(species_name)

        if merged_df is not None and not merged_df.empty:
            # Create species-specific results directory
            species_results_dir = os.path.join(RESULTS_DIR, species_name)
            os.makedirs(species_results_dir, exist_ok=True)

            # Define output file path
            output_filename = f"{species_name}_tpm_matrix.csv"
            output_filepath = os.path.join(species_results_dir, output_filename)

            try:
                merged_df.to_csv(output_filepath)
                print(f"Successfully saved TPM matrix for {species_name} to: {output_filepath}")
            except Exception as e:
                print(f"Error saving CSV for {species_name}: {e}")
        elif merged_df is not None and merged_df.empty:
             print(f"Warning: Merged DataFrame for {species_name} is empty. No CSV file will be saved.")


    print("\nScript finished.")

if __name__ == "__main__":
    main()
