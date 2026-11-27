#!/bin/bash
# Purpose: Trim paired-end `*_1.fastq`/`*_2.fastq` reads with Trimmomatic 0.39 and TruSeq3 adapters.
# Manuscript (Peleke et al., 2026): Preprocessing of paired-end RNA-seq reads (using `TruSeq3-PE.fa`)
#   before `kallisto quant` (`*_trimmed_paired_{for,rev}.fastq`).
# Usage: Run inside the directory holding raw `*_1.fastq`/`*_2.fastq` pairs after setting `ADAPTERS`.
# How to change or extend: Adjust `SLIDINGWINDOW`/`MINLEN` or adapter FASTA path (`TruSeq3-PE.fa`).
# Set adapter FASTA path (defaults to TruSeq3-PE.fa bundled alongside this script)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ADAPTERS="${ADAPTERS:-$SCRIPT_DIR/TruSeq3-PE.fa}"

# Function to perform trimming for paired-end samples
trim_sample_paired() {
  input1="$1"
  input2="$2"
  output_forward_paired="${input1%_1.fastq}_trimmed_paired_for.fastq"
  output_forward_unpaired="${input1%_1.fastq}_forward_unpaired.fastq"
  output_reverse_paired="${input2%_2.fastq}_trimmed_paired_rev.fastq"
  output_reverse_unpaired="${input2%_2.fastq}_reverse_unpaired.fastq"
  log_file="${input1%_1.fastq}_trimmomatic.log"

  trimmomatic PE -threads 8 "$input1" "$input2" \
    "$output_forward_paired" "$output_forward_unpaired" \
    "$output_reverse_paired" "$output_reverse_unpaired" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10:2:True \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > "$log_file" 2>&1
}

# Function to trim all FASTQ files in the current directory
trim_all_samples() {
  for input1 in ./*_1.fastq; do
    input2="${input1/_1.fastq/_2.fastq}"
    if [[ -f "$input1" && -f "$input2" ]]; then
      trim_sample_paired "$input1" "$input2"
    else
      echo "Skipping pair: $input1 and $input2 (one or both files not found)"
    fi
  done
}

# Ensure the adapter file exists
if [ ! -f "$ADAPTERS" ]; then
  echo "Error: Adapter file '$ADAPTERS' not found."
  exit 1
fi

# Call the function to trim all samples in the current directory
trim_all_samples

