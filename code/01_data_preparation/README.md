# Data Preparation Pipeline (`code/01_data_preparation/`)

Scripts that turn raw SRA/ENA RNA-seq runs and Helixer genome annotations into the TPM matrices, binary expression labels, threshold reports, and paralog-filtered validation gene sets committed under [`data/`](../../data/README.md) and [`results/label_thresholds_20261126/`](../../results/label_thresholds_20261126/label_threshold_old_vs_new_20261126.tsv).

> **Note on portability:** The 4 intermediate library-building scripts (`00_built_read_count_database_dPHY.0.py`, `01_built_libraries_log10Mean_dPHY.0.py`, `harmonize_libraries_dPHY.0.py`, `02_classify_libaries_dPHY.0.3.py`) default to `../data/` constants at the top of each file; adjust paths before re-running directly. The remaining scripts (`filterSRAs.py`, `downloadExps.py`, `SRA_fast-dump.v2.sh`, `trimm_PE.sh`, `process_trimmed_reads.sh`, `kallisto_SE.sh`, `run_my_diamond_jobs.sh`) use repo-relative `data/` or `EXTERNAL/` paths, while the standalone publication scripts (`downloadExps.1.py`, `make_labels_dPHYlogic_meanTPM.0.py`, `make_label_reports_dPHYlogic_meanTPM.0.py`, `label_threshold_summary.{0,1}.py`, `prepare_non_homologous_validation_sets_dPHY.8.py`) accept explicit CLI paths.

---

## 1. RNA-seq Retrieval, Trimming, QC, and Quantification

| Script | Purpose | Manuscript role | Usage / Notes |
|---|---|---|---|
| `filterSRAs.py` | Filters each species' `SraRunTable` to paired-end RNA-seq runs from root, shoot, leaf, or whole plant (`SraRunTable_filtered_data_V3.csv`). | Run selection for 11 species (note: final run set also includes 22 single-end runs listed in `data/rnaseq_runs_used.tsv`). | Edit `root_directory` and run `python3 code/01_data_preparation/filterSRAs.py`. |
| `downloadExps.py` | Original `prefetch` helper covering 4 initial species with hardcoded run tuples. | Early acquisition prototype superseded by `.1.py`. | Edit `base_dir`; run `python3 code/01_data_preparation/downloadExps.py`. |
| `SRA_fast-dump.v2.sh` | Original shell wrapper running `fastq-dump --split-3` on a text list of SRR accessions. | Paired/single-end FASTQ extraction helper. | `bash code/01_data_preparation/SRA_fast-dump.v2.sh -s <sra_list> -p <project_dir>` |
| `downloadExps.1.py` | Reads `data/rnaseq_runs_used.tsv` and runs SRA Toolkit `prefetch` + `fastq-dump --split-3` for all 154 runs across 13 species. | Complete reproducible retrieval of all RNA-seq reads entering the 2026 TPM matrices. | `python3 code/01_data_preparation/downloadExps.1.py <base_dir> [--runs data/rnaseq_runs_used.tsv] [--species ...] [--tissues ...] [--threads 4]` |
| `trimm_PE.sh` | Trimmomatic 0.39 paired-end adapter/quality trimming using `TruSeq3-PE.fa`. | Trims paired-end FASTQ reads to `*_trimmed_paired_{for,rev}.fastq`. (Single-end runs in *S. lycopersicum*, *S. bicolor* root, and *Z. mays* were trimmed to `*_trimmed.fastq` externally.) | Run inside a raw-FASTQ folder; verify `ADAPTERS` / `TruSeq3-PE.fa` path. |
| `process_trimmed_reads.sh` | Runs FastQC (`-t 8`) and MultiQC on trimmed FASTQ files. | Post-trimming read quality verification. | `bash code/01_data_preparation/process_trimmed_reads.sh -i <trimmed_dir> -o <qc_out_dir>` |
| `kallisto_SE.sh` | Builds a kallisto transcript index and runs single-end `kallisto quant --single -l 200 -s 20`. | Single-end quantification (`kallisto 0.50.1`). Paired-end calls are recorded per run in `data/rnaseq_runs_used.tsv` (`kallisto_call`). | `bash code/01_data_preparation/kallisto_SE.sh -i <trimmed_dir> -r <cds.fa> -o <out_dir>` |

---

## 2. TPM Matrices, Normalization, and Expression Labeling

| Script | Purpose | Manuscript role | Usage / Notes |
|---|---|---|---|
| `00_built_read_count_database_dPHY.0.py` | Aggregates per-run `abundance.tsv` files into `<species>_tpm_matrix.csv`. | Builds `data/<species>/<species>_tpm_matrix.csv.gz` across the 13 species. | Configure `BASE_DATA_DIR` / `RESULTS_DIR`; run inside `bin/` layout or update constants. |
| `01_built_libraries_log10Mean_dPHY.0.py` | Computes `log10(mean TPM + 1)` per tissue (`leaf`, `root`) and across all tissues. | Computes expression metric shown in **Figure 1b** and used for binary labeling. | Configure `BASE_INPUT_DIR`; reads `*_tpm_matrix.csv`. |
| `harmonize_libraries_dPHY.0.py` | Standardizes normalized CSVs to `gene_id, log10MeanTPM, target` and strips trailing isoform/version suffixes (`.1`, `.2`). | Prepares normalized per-tissue tables for `02_classify_libaries_dPHY.0.3.py`. | Configure `BASE_RESULTS_DIR`; overwrites normalized files in place. |
| `02_classify_libaries_dPHY.0.3.py` | Implements `classify_genes()` (deepCRE2.0 class rule) and `generate_report_and_plot()` (threshold histograms & text reports). | Core classification engine imported unchanged via `importlib` by `make_labels_*` and `make_label_reports_*`. | Standalone `--norm_type mean --file_stem leaf --q_threshold 0.75 --random_seed 2026`, or imported as library module. |
| `make_labels_dPHYlogic_meanTPM.0.py` | Runs `classify_genes()` (`q = 0.75`, `seed = 2026`) across all 13 species × leaf/root. | Generates `data/<species>/*_labeled_dPHYlogic.csv.gz` (`1` = high `> P75`, `0` = low matched count with `<= 50%` zero-TPM, `2` = unused) and `data/label_summary.tsv`. | `python3 code/01_data_preparation/make_labels_dPHYlogic_meanTPM.0.py <normalized_root> <out_root> [--q 0.75] [--seed 2026]` |
| `make_label_reports_dPHYlogic_meanTPM.0.py` | Verifies stored `target` labels against `classify_genes()` and runs `generate_report_and_plot()`. | Generates `data/<species>/label_reports/{reports,plots}/*` across all 26 sets. | `python3 code/01_data_preparation/make_label_reports_dPHYlogic_meanTPM.0.py data [--q 0.75] [--seed 2026]` |
| `label_threshold_summary.0.py` | Computes old (`low <= P25, high >= P75`) vs. 2026 deepCRE2.0 thresholds and TPM equivalents (`10**x - 1`) across label files. | Intermediate metrics for **Table 2**. | `python3 code/01_data_preparation/label_threshold_summary.0.py <labels_root> <summary_tsv>` |
| `label_threshold_summary.1.py` | Formats `.0.py` output into the side-by-side **Table 2** comparison table. | Produces `results/label_thresholds_20261126/label_threshold_old_vs_new_20261126.tsv`. | `python3 code/01_data_preparation/label_threshold_summary.1.py <summary_tsv_from_.0> <out_tsv>` |

---

## 3. Gene Models and Paralog-Filtered Validation Sets (`v8`)

* **Helixer** (`land_plant_v0.3_a_0080`; md5 `f0e00efcbea83c66b69258d11119a691`) called on each RefSeq `.fna` assembly produces `data/<species>/*.helixer.gff.gz`; **AGAT v1.4.0** converts them to `*.helixer.agat.gtf.gz` for 3,020 bp window extraction.

| Script | Purpose | Manuscript role | Usage / Notes |
|---|---|---|---|
| `generate_diamond_commands.py` | Generates shell commands for per-species DIAMOND `makedb` and self-`blastp`. | Created the template for `run_my_diamond_jobs.sh`. | `python3 code/01_data_preparation/generate_diamond_commands.py --proteome_dir <dir> --diamond_db_dir <db> --diamond_out_dir <out>` |
| `run_my_diamond_jobs.sh` | Runs DIAMOND 2.1.9 self-`blastp` (`--ultra-sensitive --max-hsps 1 --max-target-seqs 25 --evalue 1e-3`). | Produces `EXTERNAL/diamond/BLASTdP_<sp>_to_<sp>` all-vs-all protein alignments. | Decompress proteomes (or regenerate commands via `generate_diamond_commands.py`) and run `bash code/01_data_preparation/run_my_diamond_jobs.sh` from repo root. |
| `prepare_non_homologous_validation_sets_dPHY.7.py` | Earlier validation-gene selector treating any non-self DIAMOND hit (`E <= 1e-3`) as a paralog with tier-3 fallback fill. | Legacy `v7` baseline kept for documented comparison in `.8.py`. | Superseded by `.8.py`. |
| `prepare_non_homologous_validation_sets_dPHY.8.py` + `validation_v8/species_config_v8.json` | Selects chromosome-stratified validation genes free of inter-chromosomal paralogs (`pident >= 50%`, reciprocal coverage `>= 50%`, `E <= 1e-10`, no tier-3 fill, `seed = 2026`). | Produces `data/validation_set_dCRE2026.v8.p50c50e10.pickle` used in Methods, training, Figure 2, and Supplementary Figure 1. | `python3 code/01_data_preparation/prepare_non_homologous_validation_sets_dPHY.8.py --species_config code/01_data_preparation/validation_v8/species_config_v8.json --output <out.pickle>` |

---

## How to Change or Extend

* **Add a species:** Append its runs to `data/rnaseq_runs_used.tsv`, run steps 1–2, add its proteome and DIAMOND path to `validation_v8/species_config_v8.json`, and re-run `prepare_non_homologous_validation_sets_dPHY.8.py`.
* **Adjust stringency or class quantile:** Change `--min_pident`, `--min_cov`, `--max_evalue`, or `--n_per_species` on `.8.py`, or `--q` / `--seed` on `make_labels_dPHYlogic_meanTPM.0.py`.
