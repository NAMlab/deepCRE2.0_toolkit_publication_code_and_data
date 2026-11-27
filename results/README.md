# Evaluation and Threshold Result Tables (`results/`)

Pre-computed summary tables produced by `code/01_data_preparation/` and `code/03_evaluation/`. All figure scripts in `code/04_figures/` and workbook builders (`build_source_data_fig2.5.py`, `build_supplementary_data1.1.py`) read directly from this directory so figures and supplementary workbooks can be regenerated without re-running cluster inference.

## Directory Structure

### `results/label_thresholds_20261126/`
* **`label_threshold_old_vs_new_20261126.tsv`** — Produced by `code/01_data_preparation/label_threshold_summary.1.py` (from `.0.py` output). Side-by-side comparison per species and tissue of the Peleke et al. (2024) expression cut-offs (`low <= P25`, `high >= P75`) vs. the 2026 deepPHY-logic rule (`high > P75`, `low` matched in count with `<= 50%` zero-TPM genes). Supports **Table 2** and Methods.

### `results/evaluation/final/`
Produced by `code/03_evaluation/s2v8_validation_metrics.2.py`, `supptable1_fold_counts.0.py`, and `crossspecies_matrices.2.py`:
* `S2v8_per_fold.tsv` — Per-fold metrics for all 231 converged/final single-species models (`v8` paralog-filtered balanced draws and full held-out chromosome metrics) + `deployed` selection flag (highest validation `auROC_v8` per species × tissue). Feeds Figure 2b–c, Supplementary Figure 1b–c, and Supplementary Data 1 (`per_fold_metrics`).
* `S2v8_per_species.tsv`, `S2v8_per_tissue.tsv` — Species × tissue and tissue-level summaries of single-species validation performance.
* `deployed_models_leaf.tsv`, `deployed_models_root.tsv` — The 13 leaf and 13 root deployed single-species fold filenames (`26` total). Feeds `build_supplementary_data1.1.py` (`deployed_models` sheet) and `crossspecies_matrices.2.py`.
* `SuppTable1_per_fold_counts.tsv` — Per-fold training and validation gene counts alongside `v8` metrics (`supptable1_fold_counts.0.py`) for **Supplementary Table 1** and Supplementary Data 1 (`fold_sizes` sheet).
* `crossspecies_mean_{leaf,root}_class{0,1}.tsv` — 13 × 13 matrices of mean predicted probability of high expression across all deployed model species (rows) × target species (columns), separated by true expression class (`0` = low, `1` = high). Feeds Figure 2d–e and Supplementary Figure 1d–e.

### `results/evaluation/msr/`
Produced by `code/03_evaluation/msr_eval_dPm.2.py` and `msr_matrices.0.py`:
* `msr_eval_{leaf,root}_full.tsv`, `msr_eval_{leaf,root}_summary.tsv` — Leave-one-species-out evaluation of the 52 multi-species (`all`, `dicots`, `monocots`) models across held-out, seen, and out-of-clade target species (`100` balanced draws per comparison). Feeds Figure 2f–s, Supplementary Figure 1f–s, and Supplementary Data 1 (`msr_evaluation`, `msr_summary`).
* `suppfig_msr_mean_{leaf,root}_{all,dicots,monocots}_class{0,1}.tsv` — 12 heatmap matrices (all-species, dicots, monocots × leaf/root × true class `0`/`1`) plotted in Figure 2i/m/n/r/s and Supplementary Figure 1i/m/n/r/s.

### `results/evaluation/ssc_control/`
Produced by `code/03_evaluation/ssc_control_table.0.py`:
* `ssc_control_per_fold.tsv`, `ssc_control_per_species.tsv` — Shuffled-sequence control vs. native-sequence validation metrics across all 231 chromosome folds (`13` species × `2` tissues, including *H. vulgare*). Feeds Supplementary Data 1 (`shuffled_control` sheet) and Methods.

### `results/evaluation/tissue_cross/`
Produced by `code/03_evaluation/tissue_*.py`:
* `tissue_cross_full.tsv`, `tissue_cross_own_species.tsv` — Same-tissue vs. cross-tissue evaluation of all 26 deployed SSR models restricted to genes scored in both tissues (`tissue_cross_eval.1.py`).
* `tissue_cross_msr_full.tsv`, `tissue_cross_msr_switch_direction.tsv` — Same-tissue vs. cross-tissue evaluation and switcher contrasts for MSR models (`tissue_cross_msr.0.py`).
* `tissue_label_composition.tsv` — Per-species counts of genes stable, swapped, or tissue-specific/intermediate between leaf and root (`tissue_label_composition.1.py`).
* `tissue_switchers.tsv`, `tissue_switch_direction.tsv`, `tissue_switch_performance.tsv`, `tissue_switch_permutation.tsv` — Directional probability contrast `P(leaf model) - P(root model)`, balanced classification metrics, and permutation tests on genes whose expression class flips between leaf and root (`tissue_cross_switchers.1.py`, `tissue_switch_direction.1.py`, `tissue_switch_performance.0.py`). Feeds **Supplementary Figure 2**.

## How to Regenerate or Extend

Pass `results/evaluation` directly to downstream figure or workbook scripts, or point `code/03_evaluation/` scripts at a new `EXTERNAL/predictions/` run and write outputs to a parallel results folder (e.g. `results/evaluation_custom/`) before passing that folder to `S2_dir` / `build_supplementary_data1.1.py` / `build_source_data_fig2.5.py`.
