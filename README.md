[![DOI](https://zenodo.org/badge/github/1474269400.svg)](https://doi.org/10.5281/zenodo.24597674)
# deepCRE Toolkit (Peleke et al., 2026): Publication Code, Data, Models, and Results

Repository accompanying Peleke et al. (2026) (*"dAccessible multi-species prediction of plant gene cis-regulation with deepCRE"*). It collects the curated expression data, gene annotations, training configurations, trained model archives, evaluation tables, and figure/note scripts behind the deepCRE models served at https://deepcre.ipk-gatersleben.de (13 flowering plant species × leaf and root; 232 single-species [SSR], 217 archived shuffled-sequence control [SSC], and 52 multi-species [MSR] models, plus the 45 Peleke et al. 2024 reference models).

All Code, Data and Models have been reviewed and are authored by Simon M. Zumkeller.

---

## Repository Layout

| Path | Contents | Documentation |
|---|---|---|
| `data/` | 13 species folders (Helixer GFF3 and AGAT GTF gene models, Ensembl Plants release 59 GFF3s, kallisto TPM matrices, deepCRE2.0 binary labels, class-threshold reports/histograms, Helixer proteomes, SRA metadata), run manifests (`genomes.tsv`, `rnaseq_runs_used.tsv`, `rnaseq_runs.tsv`), and paralog-filtered validation gene pickles. | [`data/README.md`](data/README.md) |
| `code/01_data_preparation/` | Upstream RNA-seq download, Trimmomatic trimming, FastQC/MultiQC, kallisto quantification, TPM matrix construction, log10(mean TPM + 1) normalization, deepPHY binary labeling & histograms, Table 2 threshold comparison, and DIAMOND paralog-filtered validation set generation (`v8`). | [`code/01_data_preparation/README.md`](code/01_data_preparation/README.md) |
| `code/02_training/` | Seed-and-restart wrapper (`train_models_seed_restart.{1,2}.py`), SGE cluster job scripts and JSON configs for SSR, MSR (all-species, dicots, monocots), SSC controls, model selection/ID assignment, and the *A. thaliana* TAIR10 annotation-frame exception. | [`code/02_training/README.md`](code/02_training/README.md) |
| `code/03_evaluation/` | Configs and cluster scripts for held-out-chromosome and cross-species inference, plus balanced validation metrics, Supplementary Table 1 fold counts, Supplementary Data 1 workbook builder, web-model manifest builder, MSR evaluation, and cross-tissue switcher analyses. | [`code/03_evaluation/README.md`](code/03_evaluation/README.md) |
| `code/04_figures/` | Publication figure scripts for Figure 1, Figure 2 (leaf), Supplementary Figure 1 (root), Supplementary Figure 2 (cross-tissue), and the Source Data Excel workbook. | [`code/04_figures/README.md`](code/04_figures/README.md) |
| `code/05_supplementary_figures/` | TensorFlow 2.10 inference/DeepSHAP scripts and pre-computed tables reproducing Supplementary Figures 5 (promoter swaps + ZmUBI1 + conserved-block mutations) and 6 (RAP2.12 haplotype/permutation analyses across 2024 and 2026 models). | [`code/05_supplementary_figures/README.md`](code/05_supplementary_figures/README.md) |
| `code/06_annotation_frame/` | Scripts and summary TSVs for Supplementary Note 4 (quantifying 5' TSS frame differences and prediction accuracy between Helixer training windows and TAIR10/RefSeq curated annotations across species). | [`code/06_annotation_frame/README.md`](code/06_annotation_frame/README.md) |
| `code/06_benchmarking/` | Benchmark of deepCRE against PhytoExpr and PlantCaduceus (PlantCAD2) on the *A. thaliana* and *Z. mays* genomes (Figure 3 and Supplementary Figure 3): CWL workflows, Docker/Singularity container definitions, environments and setup scripts; the model code is included as git submodules in `workflows/` (`deepcre`, `phytoexpr`, `plantcad`, `ModelPipelineIO`, `PanGenomeTools`). | [`code/06_benchmarking/README.md`](code/06_benchmarking/README.md) |
| `results/` | Final evaluation tables (`final/`, `msr/`, `ssc_control/`, `tissue_cross/`) and Table 2 expression-threshold comparison (`label_thresholds_20261126/`). | [`results/README.md`](results/README.md) |
| `models/` | Trained `.h5` models in `ssr/`, `ssc/`, `msr/`, `peleke2024/`; `model_id_mapping.tsv` (515 2026 models), and `peleke2024_models.tsv` (45 models from Peleke et al., 2024). | [`models/README.md`](models/README.md) |
| `tutorial/` | Web-tool tutorial inputs (`File_1`–`File_3`, corresponding to Supplementary files 1–3) and tutorial result CSVs (`Tutorial_table{1,2,3}_*.csv`). | [`tutorial/README.md`](tutorial/README.md) |
| `code/06_benchmarking/workflows/deepcre/` | Git submodule (shared with the benchmark) pointing to https://github.com/NAMlab/deepCRE2.0 (formerly `deepCRE_reimplemented`) at commit `a37db45`, containing the unmodified core CNN training and prediction code (`train_models.py`, `deepcre_crosspredict.py`). The published models were trained at `19f6db7`, two commits earlier; `a37db45` only skips genes whose flanks reach past a contig end (instead of stopping with an error) and adds an error check to `deepcre_interpret.py`. | — |

---

## Benchmark (`code/06_benchmarking/`)

The benchmark of deepCRE against PhytoExpr and PlantCaduceus (PlantCAD2) for Figure 3 and Supplementary Figure 3. Installation and run instructions are in [`code/06_benchmarking/README.md`](code/06_benchmarking/README.md).

**Included in this repository**

| Path in `code/06_benchmarking/` | Contents |
|---|---|
| `README.md` | Requirements, installation (cwltool, Apptainer/Singularity or Docker, NVIDIA container toolkit) and how to run the workflows. |
| `environment.yml` | Conda environment with `cwltool`. |
| `setup-apptainer.sh`, `setup-docker.sh` | Build the model containers (Apptainer SIF files or Docker images). |
| `setup-run-aliases.sh` | Shell aliases (`mpl-run`, with Docker/Apptainer and parallel options) for running a workflow. |
| `workflows/cwl/deepcre/` | CWL tool `predict_crosspredict.cwl` (runs `deepcre_crosspredict.py` from the `workflows/deepcre` submodule), `Dockerfile`, Apptainer `Container.def`, `environment.yml`. |
| `workflows/cwl/phytoexpr/` | CWL tool `predict.cwl` with `predict.py`, `Dockerfile`, `Container.def`, `environment.yml`. |
| `workflows/cwl/plantcad/` | CWL tools `tokenize.cwl` and `predict_loraexpression.cwl`, `Dockerfile`, `Container.def`, `environment.yml`, `additional-requirements.txt`. |
| `workflows/cwl/deepcis/` | CWL tool `predict.cwl` with `predict.py`, `Dockerfile`, `Container.def`, `environment.yml`. |
| `workflows/deepcre` | Git submodule: https://github.com/NAMlab/deepCRE2.0 at `a37db45` (also used by `code/02_training/` and `code/03_evaluation/`). |
| `workflows/phytoexpr` | Git submodule: https://github.com/NAMlab/PhytoExpr-inofficial at `5fda808`. |
| `workflows/plantcad` | Git submodule: https://github.com/plantcad/plantcad at `bc2503d`. |
| `workflows/PanGenomeTools` | Git submodule: https://github.com/pfennigt/PanGenomeTools at `179c7c0`. |
| `workflows/ModelPipelineIO` | Git submodule: https://github.com/pfennigt/ModelPipelineIO at `55b64b6` (repository not public; access required). |

**Not included in this repository:** the `runs/` folders that `code/06_benchmarking/README.md` refers to (per-model `run.cwl`/`run.yml`, the predictions shown in the manuscript in each run's `dataset/` folder, `logs/`) and the analysis notebook that draws the figures (`runs/analysis-deepcreComparison/analysis.ipynb`). They are part of the ARC https://git.nfdi4plants.org/szymanskilab/ModelPipeline_DeepCRETK (access required).

---

## General Conventions

* **Run from repository root:** All analysis and figure scripts assume relative paths starting from the root of this repository (unless explicitly noted as requiring staging alongside the web-app `lib/` modules in `code/05_supplementary_figures/` and `code/06_annotation_frame/`).
* **External inputs (`EXTERNAL/`):** Large non-git artifacts (prediction CSVs from cluster inference, raw cluster log directories, RefSeq GFFs for annotation-frame checks, DIAMOND tabular outputs, web-model manifest, and output directories) use explicit `EXTERNAL/<category>/...` placeholders annotated with `# EXTERNAL` comments and documented in [`EXTERNAL_INPUTS.md`](EXTERNAL_INPUTS.md).
* **Script versioning:** Scripts carry numeric suffixes (`name.N.py`). Only the version used for the final manuscript (plus any earlier version imported at runtime or compared in documentation) is included.

---

## Manuscript Item → Code → Data Map

| Manuscript item | Code | Primary inputs → outputs |
|---|---|---|
| Expression labeling & Table 2 (2024 vs 2026 low-class definition) | `code/01_data_preparation/make_labels_dPHYlogic_meanTPM.0.py`, `make_label_reports_dPHYlogic_meanTPM.0.py`, `label_threshold_summary.{0,1}.py` | `data/<species>/*_tpm_matrix.csv.gz` → `data/<species>/*_labeled_dPHYlogic.csv.gz`, `data/label_summary.tsv`, `results/label_thresholds_20261126/` |
| Paralog-filtered validation sets (Methods) | `code/01_data_preparation/prepare_non_homologous_validation_sets_dPHY.8.py` + `validation_v8/species_config_v8.json` | Helixer proteomes (`data/<species>/`) + `EXTERNAL/diamond/` → `data/validation_set_dCRE2026.v8.p50c50e10*.pickle` |
| Model training (SSR, MSR, SSC, TAIR10 exception) | `code/02_training/train_models_seed_restart.{1,2}.py`, `s2v8_select_models.1.py`, `assign_model_ids.0.py`, job scripts & configs | `data/` + RefSeq `.fna` → `models/{ssr,ssc,msr}/` + `code/02_training/final_models.tsv` + `models/model_id_mapping.tsv` |
| Held-out, cross-species, MSR, SSC & cross-tissue evaluation; Supplementary Table 1; Supplementary Data 1 | `code/03_evaluation/` (`s2v8_validation_metrics.2.py`, `crossspecies_matrices.2.py`, `msr_eval_dPm.2.py`, `msr_matrices.0.py`, `ssc_control_table.0.py`, `tissue_*.py`, `supptable1_fold_counts.0.py`, `build_supplementary_data1.1.py`) | `EXTERNAL/predictions/`, `data/` → `results/evaluation/` & Supplementary Data 1 Excel workbook |
| Figure 1 (web tool workflow) | `code/04_figures/build_updated_fig1.py` | `EXTERNAL/figure_inputs/` → `EXTERNAL/figure_output/figures/` |
| Figure 2 (leaf), Supplementary Figure 1 (root), Supplementary Figure 2 (cross-tissue), Source Data | `code/04_figures/make_fig2_comprehensive.{0,2}.py`, `make_suppfig1_root.2.py`, `make_suppfig2_tissue.7.py`, `make_pub_figures_split.7.py`, `build_source_data_fig2.5.py` | `results/evaluation/` → `EXTERNAL/figure_output/` (or `$FIG_OUT_DIR`) + Source Data `.xlsx` |
| Supplementary Figures 5 (promoter swaps) and 6 (RAP2.12 variants) | `code/05_supplementary_figures/` (11 inference/DeepSHAP scripts + `make_supp_figures_full.{0,15}.py`, `make_suppfig6_classic.{0,14}.py`) | `tutorial/`, `models/`, TSV/NPZ intermediates in `code/05_supplementary_figures/` → figure PNGs + legend `.txt` |
| Supplementary Note 4 (annotation frame comparison) | `code/06_annotation_frame/annotation_shift_run{A,B,C,D}.0.py`, `annotation_frame_runE.1.py`, `annotation_frame_runF.0.py` | `data/`, `models/`, `EXTERNAL/refseq_gff/` → `code/06_annotation_frame/*.tsv` |
| Figure 3 and Supplementary Figure 3 (benchmark against PhytoExpr and PlantCAD2) | `code/06_benchmarking/` (CWL workflows and model submodules) | *A. thaliana* and *Z. mays* genomes and expression data → per-model predictions and benchmark metrics; see [`code/06_benchmarking/README.md`](code/06_benchmarking/README.md) |
| *A. thaliana* curated-window TAIR10 models | `code/02_training/tair10_exception/` | `arabidopsis_thaliana_{leaf,root}_*_TAIR10.csv`, `validation_set_*.tair10.pickle` → `Atha_TAIR10_*` models |

> **Reproducibility check:** Running `make_fig2_comprehensive.2.py`, `make_suppfig1_root.2.py`, and `make_suppfig2_tissue.7.py` from `results/` reproduces the plotted values of Figure 2, Supplementary Figure 1, and Supplementary Figure 2 byte for byte, and `build_source_data_fig2.5.py` and `build_supplementary_data1.1.py` reproduce every sheet of the Source Data and Supplementary Data 1 workbooks.

---

## Environment Setup & Quick Start

1. **Clone the `deepCRE2.0` submodule (only needed if re-training or running `deepcre_crosspredict.py`):**
   ```bash
   git submodule update --init code/06_benchmarking/workflows/deepcre
   git -C code/06_benchmarking/workflows/deepcre checkout 19f6db7   # optional: the exact commit the published models were trained with
   cp code/02_training/train_models_seed_restart.*.py code/06_benchmarking/workflows/deepcre/src/deepCRE/
   ```
2. **Decompress annotations/labels/matrices and download RefSeq genome FASTAs (only needed for training/inference):**
   ```bash
   gunzip -k data/*/*.gz
   bash data/download_genomes.sh
   ```
3. **Trained HDF5 models** are in `models/ssr/`, `models/ssc/`, `models/msr/` and `models/peleke2024/` (see `models/README.md`).
4. **Install the analysis Python environment (Python 3.11) to re-run evaluation, figure, and workbook scripts:**
   ```bash
   pip install -r requirements_analysis.txt
   ```
   Example (generate Figure 2, Supplementary Figure 1, Supplementary Figure 2, and workbooks into `EXTERNAL/figure_output/pub_out` without editing code):
   ```bash
   FIG_OUT_DIR=EXTERNAL/figure_output/pub_out python3 code/04_figures/make_fig2_comprehensive.2.py
   FIG_OUT_DIR=EXTERNAL/figure_output/pub_out python3 code/04_figures/make_suppfig1_root.2.py
   FIG_OUT_DIR=EXTERNAL/figure_output/pub_out python3 code/04_figures/make_suppfig2_tissue.7.py
   mkdir -p EXTERNAL/figure_output/pub_out
   python3 code/03_evaluation/build_supplementary_data1.1.py \
     results/evaluation results/evaluation/ssc_control/ssc_control_per_fold.tsv EXTERNAL/figure_output/pub_out/Supplementary_Data_1.xlsx
   ```

Training (`deepCRE2.0`) and TensorFlow inference/DeepSHAP scripts (`code/05_supplementary_figures/`, `code/06_annotation_frame/`) require Python 3.8 + TensorFlow 2.10 (see `code/06_benchmarking/workflows/deepcre/pyproject.toml`). Scripts importing `from lib.utils import ...` additionally expect the `lib/` package from https://github.com/NAMlab/streamlit_deepcre_app.

---

## Using This Repository as a Lab Template

This repository is structured so new lab members can adapt and extend individual stages without rebuilding the whole pipeline:

1. **Add a new species or tissue:**
   * Add the species and RefSeq assembly accession to `data/genomes.tsv` and run `bash data/download_genomes.sh`.
   * Prepare RNA-seq TPM matrices and binary labels following [`code/01_data_preparation/README.md`](code/01_data_preparation/README.md) (`make_labels_dPHYlogic_meanTPM.0.py`), and add the new species entry to `code/01_data_preparation/validation_v8/species_config_v8.json` before running `prepare_non_homologous_validation_sets_dPHY.8.py`.
   * Create a training JSON modeled on `code/02_training/ssr/configs/01_Atha_leaf.json` and train with `train_models_seed_restart.2.py`.
   * To include the new species in evaluation matrices and figures, add its 4-letter code, label, clade topology, and chromosome ordering to `ORDER` / `LABELS` / `TOPO_ALL` / `CLADES` in `code/03_evaluation/` and `code/04_figures/make_pub_figures_split.7.py`.
2. **Define a new multi-species (MSR) training subset:**
   * Add a JSON config in `code/02_training/msr/configs/` (copying e.g. `01_MSR_all_leaf.json` and editing the `"species"` list), train with `train_models_seed_restart.2.py`, generate cross-predictions via `code/03_evaluation/crosspred_msr/`, and evaluate with `msr_eval_dPm.2.py` and `msr_matrices.0.py`.
3. **Modify expression class thresholds or quantile rules:**
   * Pass custom `--q` / `--seed` flags to `make_labels_dPHYlogic_meanTPM.0.py` or adapt `classify_genes()` in `02_classify_libaries_dPHY.0.3.py`, then inspect shifts via `label_threshold_summary.0.py` and `label_threshold_summary.1.py`.
4. **Benchmark models on an alternative gene annotation frame:**
   * Set `RUNE_SPECIES` and `RUNE_REFSEQ_DIR` (see `code/06_annotation_frame/annotation_frame_runE.1.py`) or point a config JSON to an alternative GTF and pickle key (as in `code/02_training/tair10_exception/`).
5. **Re-style or extend downstream figures and workbooks:**
   * Redirect figure outputs via `FIG_OUT_DIR=/path/to/out` or pass CLI output directories to `make_supp_figures_full.15.py`, `make_suppfig6_classic.14.py`, `build_source_data_fig2.5.py`, and `build_supplementary_data1.1.py`.

---

## AI-Assisted Editing and Authorization

The deepCRE pipeline, models, and analyses were developed by the study authors. AI coding assistants (Claude Code Opus5 and Opus5.5 and Google Antigravity) were used selectively to help refactor helper/evaluation scripts, adapt configurations, and structure repository documentation; the codebase and the data is not AI-generated. All code, configurations, data tables, and documentation in this repository were verified against the manuscript outputs and reviewed and authorized by **S.M.Z.**
