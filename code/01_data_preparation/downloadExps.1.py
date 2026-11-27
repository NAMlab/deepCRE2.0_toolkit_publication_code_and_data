"""
Purpose:
    Download SRA/ENA RNA-seq runs and convert them to FASTQ using SRA Toolkit
    (`prefetch` and `fastq-dump --split-3`) for all 154 runs in `data/rnaseq_runs_used.tsv`.

Manuscript (Peleke et al., 2026):
    Expanded, reproducible replacement for `downloadExps.py` + `SRA_fast-dump.v2.sh`
    covering all 13 species and leaf/root/shoot runs underlying the 2026 TPM matrices.

Usage & I/O (from repo root; requires `prefetch` and `fastq-dump` on PATH):
    python3 code/01_data_preparation/downloadExps.1.py <base_dir> \
        [--runs data/rnaseq_runs_used.tsv] [--species zea_mays,...] [--tissues leaf,root] [--threads 4]
    - Outputs under `<base_dir>/<species>/data/<tissue>/raw/`: `.sra`, `_1.fastq`/`_2.fastq`
      (paired-end), or `<RUN>.fastq` (single-end, matching `kallisto_call` in `rnaseq_runs_used.tsv`).

How to change or extend:
    Filter subsets via `--species` / `--tissues`, increase `--threads`, or append new runs
    to a copy of `data/rnaseq_runs_used.tsv` and pass `--runs`.
"""
import argparse
import concurrent.futures
import csv
import os
import subprocess


def run_command(command):
    """Run a shell command and return its output."""
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return result.stdout, result.stderr, result.returncode


def download_run(base_dir, species, tissue, run):
    """Download a sequencing run and convert it to FASTQ in the appropriate directory."""
    raw_dir = os.path.join(base_dir, species, "data", tissue, "raw")
    os.makedirs(raw_dir, exist_ok=True)
    for command in (f"prefetch {run} -O {raw_dir}",
                    f"fastq-dump --split-3 --outdir {raw_dir} {os.path.join(raw_dir, run, run + '.sra')}"):
        stdout, stderr, returncode = run_command(command)
        if returncode != 0:
            print(f"Failed: {command}\n{stderr}")
            return False
    print(f"Downloaded {run} ({species}, {tissue}) to {raw_dir}")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base_dir")
    ap.add_argument("--runs", default="data/rnaseq_runs_used.tsv")
    ap.add_argument("--species", default=None)
    ap.add_argument("--tissues", default=None)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    with open(a.runs) as fh:
        dataset = [(r["species"], r["tissue"], r["run_accession"]) for r in csv.DictReader(fh, delimiter="\t")]
    if a.species:
        dataset = [d for d in dataset if d[0] in a.species.split(",")]
    if a.tissues:
        dataset = [d for d in dataset if d[1] in a.tissues.split(",")]
    print(f"{len(dataset)} runs")
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.threads) as executor:
        futures = [executor.submit(download_run, a.base_dir, s, t, r) for s, t, r in dataset]
        ok = sum(f.result() for f in concurrent.futures.as_completed(futures))
    print(f"{ok}/{len(dataset)} runs downloaded")


if __name__ == "__main__":
    main()
