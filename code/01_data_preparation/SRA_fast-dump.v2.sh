#!/bin/bash
# Purpose: Convert a list of SRA accessions into paired/single FASTQ files via `fastq-dump --split-3`.
# Manuscript (Peleke et al., 2026): Original FASTQ conversion wrapper paired with `downloadExps.py`
#   (superseded for all 154 runs by `downloadExps.1.py`).
# Usage: bash code/01_data_preparation/SRA_fast-dump.v2.sh -s <sra_list_file> -p <output_dir>
# How to change or extend: Swap `fastq-dump` for `fasterq-dump` inside the loop for faster multi-threaded I/O.

# Function to display usage information
usage() {
  echo "Usage: $0 -s <SRA_list_file> -p <project_name>"
  echo "Options:"
  echo "  -s <SRA_list_file>: Specify the SRA list file."
  echo "  -p <project_name>: Specify the project name."
  exit 1
}

# Default values
SRA_LIST=""
PROJECT=""

# Parse command-line arguments
while getopts "s:p:" opt; do
  case $opt in
    s)
      SRA_LIST="$OPTARG"
      ;;
    p)
      PROJECT="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      ;;
  esac
done

# Check if SRA list and project name are provided
if [ -z "$SRA_LIST" ] || [ -z "$PROJECT" ]; then
  echo "Error: Both -s (SRA list file) and -p (project name) must be specified."
  usage
fi

# Check if SRA list file exists
if [ ! -e "$SRA_LIST" ]; then
  echo "Error: SRA list file '$SRA_LIST' does not exist."
  exit 1
fi

# Create output directory if it doesn't exist
if [ ! -d "$PROJECT" ]; then
  mkdir "$PROJECT"
fi

# Process SRA accession numbers
while IFS= read -r SRR; do
  # Trim leading and trailing whitespace from each line
  SRR=$(echo "$SRR" | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//')

  # Skip empty lines and lines that start with "#" (comments)
  if [[ -n "$SRR" && "$SRR" != \#* ]]; then
    fastq-dump --split-3 "$SRR"
    if [ $? -ne 0 ]; then
      echo "Error: fastq-dump failed for accession number '$SRR'."
      continue
    fi
  fi
done < "$SRA_LIST"

# Move downloaded files to the project directory
mv SRR* "$PROJECT/"

echo "Download and organization of data for project $PROJECT is complete."

