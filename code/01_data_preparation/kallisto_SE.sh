#!/bin/bash
# Purpose: Build a kallisto transcriptome index and quantify single-end `*_trimmed.fastq` libraries.
# Manuscript (Peleke et al., 2026): Quantified single-end RNA-seq runs (`--single -l 200 -s 20`; see
#   `data/rnaseq_runs_used.tsv` for exact calls across all 154 runs).
# Usage: bash code/01_data_preparation/kallisto_SE.sh -i <trimmed_fastq_dir> -r <transcripts.fa> -o <out_dir>
# How to change or extend: Adjust `-l` / `-s` fragment length parameters or switch to paired-end arguments.

# Function to display usage information
usage() {
  echo "Usage: $0 -i <input_dir> -r <reference_fa>"
  echo "  -i <input_dir>: Specify the input directory containing the trimmed fastq files."
  echo "  -r <reference_fa>: Specify the reference FASTA file for kallisto."
  echo "  -o <output_dir>: Specify the output directory."
  exit 1
}

# Initialize variables
input_dir=""
reference_fa=""
index_file=""
output_dir=""

# Parse command line options
while getopts ":i:r:o:" opt; do
  case $opt in
    i)
      input_dir="$OPTARG"
      ;;
    r)
      reference_fa="$OPTARG"
      ;;
    o)
      output_dir="$OPTARG" 
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

# Check if input directory and reference FASTA are provided
if [ -z "$input_dir" ] || [ -z "$reference_fa" ]; then
  echo "Error: Input directory and reference FASTA are not specified."
  usage
fi

# Generate the index_file name based on the reference_fa file
reference_base=$(basename "$reference_fa")
index_file="${reference_base%.*}.index"
mkdir "${output_dir}"
# Step 1: Index the reference transcriptome
kallisto index -i "$index_file" "$reference_fa"

# Step 2: Quantify TPMs for matching pairs of files
for file in "$input_dir"/*_trimmed.fastq; do
  base=$(basename "$file")
  sample_name="${base%_trimmed.fastq}"
  kallisto quant -i "$index_file" -o "$output_dir/$sample_name" --single -l 200 -s 20 "$file"
done
