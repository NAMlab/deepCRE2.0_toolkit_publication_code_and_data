#!/usr/bin/env bash
# Purpose: Download the RefSeq genome assemblies listed in data/genomes.tsv using
#   NCBI datasets (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/).
# Manuscript (Peleke et al., 2026): Retrieves the 13 plant genome FASTA assemblies on
#   which Helixer gene models were called and from which 3,020 bp windows were extracted.
# Usage: bash data/download_genomes.sh   (runnable from repo root or any working directory)
#   Outputs: data/<species>/genome/<accession>.fna
# How to change or extend: Add a row to data/genomes.tsv (species, accession, assembly,
#   annotation, model) and re-run; existing non-empty .fna files are skipped.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
tail -n +2 "$SCRIPT_DIR/genomes.tsv" | while IFS=$'\t' read -r sp acc asm ann model; do
  out_dir="$SCRIPT_DIR/$sp/genome"
  out_fna="$out_dir/$acc.fna"
  if [ -s "$out_fna" ]; then
    echo "skip existing: $out_fna"
    continue
  fi
  mkdir -p "$out_dir"
  datasets download genome accession "$acc" --include genome --filename "$out_dir/$acc.zip"
  unzip -p "$out_dir/$acc.zip" '*.fna' > "$out_fna" && rm "$out_dir/$acc.zip"
done
