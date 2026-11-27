"""
Purpose:
    Filter per-species NCBI `SraRunTable` CSV metadata to paired-end RNA-seq libraries
    annotated for leaf, root, shoot, or whole-plant tissues (`SraRunTable_filtered_data_V3.csv`).

Manuscript (Peleke et al., 2026):
    Produced the per-species filtered run tables stored in `data/<species>/` (note:
    22 additional single-end runs present in `data/rnaseq_runs_used.tsv` were added separately).

Usage & I/O:
    python code/01_data_preparation/filterSRAs.py
    - Reads `<root_directory>/<species>/SraRunTable*.csv`
    - Writes `SraRunTable_filtered_data_V3.csv` and summary stats in each species folder.

How to change or extend:
    Set `root_directory` to your metadata folder and customize the tissue/layout filter
    conditions inside `filter_species_data()`.
"""
import pandas as pd
import os
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Define the root directory where species folders are located
root_directory = "data"

# Define a list to collect summary information for each species
summary_data = []

# Define custom logging functions
def log_processing_start(file_path):
    logging.info(f"Processing file: {file_path}")

def log_processing_end(output_file, num_rows_original, num_rows_omitted):
    num_rows_filtered = num_rows_original - num_rows_omitted
    percent_omitted = (num_rows_omitted / num_rows_original) * 100

    logging.info(f"Filtered data saved to: {output_file}")
    logging.info(f"Processing completed:")
    logging.info(f"  Original rows: {num_rows_original}")
    logging.info(f"  Filtered rows: {num_rows_filtered}")
    logging.info(f"  {percent_omitted:.2f}% rows omitted")

    # Get the species folder name
    species_folder = os.path.basename(os.path.dirname(output_file))
    
    # Append summary information for the species to the list
    summary_data.append({
        'Species': species_folder,
        'Original Rows': num_rows_original,
        'Filtered Rows': num_rows_filtered,
        'Omitted Rows': num_rows_omitted
    })

# Define the filtering logic
def filter_species_data(input_file):
    try:
        # Log processing start
        log_processing_start(input_file)

        # Read the CSV file with low_memory=False to suppress DtypeWarning
        df = pd.read_csv(input_file, low_memory=False)

        # Convert all columns to lowercase for case-insensitive matching
        df = df.apply(lambda x: x.astype(str).str.lower())

        # Define the filter conditions
        string_conditions = df.apply(lambda col: col.str.contains('rnaseq|rna|transcript*|transcript|polyA', case=False))
        paired_condition = df.apply(lambda col: col.str.contains('paired', case=False))
        tissue_condition = df.apply(lambda col: col.str.contains('root|shoot|leaf|whole', case=False))

        # Combine conditions across all columns
        overall_condition = string_conditions.any(axis=1) & paired_condition.any(axis=1) & tissue_condition.any(axis=1)

        # Count omitted rows
        num_rows_original = df.shape[0]
        num_rows_filtered = overall_condition.sum()
        num_rows_omitted = num_rows_original - num_rows_filtered

        # Log processing end and statistics
        output_file = input_file.replace('SraRunTable.txt', 'SraRunTable_filtered_data_V3.csv')
        log_processing_end(output_file, num_rows_original, num_rows_omitted)

        # Apply the filter and keep rows that satisfy the overall condition
        filtered_df = df[overall_condition]

        # Save the filtered dataframe to the new CSV file
        filtered_df.to_csv(output_file, index=False)

    except Exception as e:
        logging.error(f"Error processing file {input_file}: {e}")

# Traverse through each folder and subfolder
for dirpath, dirnames, filenames in os.walk(root_directory):
    # Check if the current folder contains SraRunTable.txt
    if 'SraRunTable.txt' in filenames:
        # Construct the full path to the input file
        input_file = os.path.join(dirpath, 'SraRunTable.txt')

        # Filter the data
        filter_species_data(input_file)

# Save the summary data to a CSV file
summary_file = os.path.join(root_directory, 'processing_summary.csv')
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(summary_file, index=False)

logging.info(f"Summary data saved to: {summary_file}")
