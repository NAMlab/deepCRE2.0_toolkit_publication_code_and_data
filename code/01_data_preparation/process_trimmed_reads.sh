#!/bin/bash
# Purpose: Run FastQC and aggregate MultiQC reports across a folder of trimmed FASTQ files.
# Manuscript (Peleke et al., 2026): Read-quality verification after adapter trimming (`trimm_PE.sh`).
# Usage: bash code/01_data_preparation/process_trimmed_reads.sh -i <trimmed_dir> -o <qc_output_dir>
# How to change or extend: Adjust FastQC thread count `-t 8` or output subdirs below.

# Usage function to display help
usage() {
    echo "Usage: $0 -i <input_directory> -o <output_directory>"
    echo "  -i <input_directory>: Directory containing trimmed FASTQ files."
    echo "  -o <output_directory>: Directory to store FastQC and MultiQC output."
    exit 1
}

# Check for the correct number of arguments
if [ "$#" -ne 4 ]; then
    usage
fi

# Parse command-line arguments
while getopts ":i:o:" opt; do
    case $opt in
        i) input_directory="$OPTARG"
        ;;
        o) output_directory="$OPTARG"
        ;;
        \?) echo "Invalid option -$OPTARG" >&2
            usage
        ;;
        :)  echo "Option -$OPTARG requires an argument." >&2
            usage
        ;;
    esac
done

# Verify input directory exists
if [ ! -d "$input_directory" ]; then
    echo "Error: Input directory does not exist."
    exit 1
fi

# Create output directory if it does not exist
if [ ! -d "$output_directory" ]; then
    mkdir -p "$output_directory"
fi

# Define FastQC output path
fastqc_output="${output_directory}/fastqc_results"
mkdir -p "$fastqc_output"

# Run FastQC on all trimmed FASTQ files
echo "Running FastQC..."
fastqc -t 8 -o "$fastqc_output" "${input_directory}"/*.fastq

# Check if FastQC ran successfully
if [ "$?" -ne 0 ]; then
    echo "FastQC did not complete successfully. Check output for errors."
    exit 1
fi

# Run MultiQC to aggregate FastQC reports
echo "Running MultiQC..."
multiqc -o "${output_directory}/multiqc_report" "$fastqc_output"

# Check if MultiQC ran successfully
if [ "$?" -ne 0 ]; then
    echo "MultiQC did not complete successfully. Check output for errors."
    exit 1
fi

echo "Quality control processing complete. Results are in: ${output_directory}/multiqc_report"

