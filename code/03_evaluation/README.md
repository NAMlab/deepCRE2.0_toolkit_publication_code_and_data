# Model Evaluation, Cross-Prediction, and Workbook Builders (`code/03_evaluation/`)

Scripts and prediction configs that run held-out and cross-species inference (`deepcre_crosspredict.py` in `deepCRE_reimplemented`), compute all balanced/unbalanced metrics under [`results/evaluation/`](../../results/README.md), build Supplementary Data 1, and generate the deployed web-tool model manifest.

---

## 1. Inference Configs & Cluster Runners

| Path | Purpose | Manuscript role | Usage / Notes |
|---|---|---|---|
| `s2v8_make_predict_config.1.py` | Generates `own_chromosome_configs/own_chromosome_predict.json` from `code/02_training/final_models.tsv` (`EXTERNAL/data_root/`). | Builds per-fold own-chromosome inference config for all 231 final SSR models. | `python3 code/03_evaluation/s2v8_make_predict_config.1.py code/02_training/final_models.tsv code/03_evaluation/own_chromosome_configs/own_chromosome_predict.json <pred_out_dir>` |
| `run_s2v8_predict.0.sh` | Runs `deepcre_crosspredict.py` locally with `$REPO/EXTERNAL/environment/python` on a config JSON in this folder. | Executes own-chromosome `S2v8` predictions (`EXTERNAL/predictions/`). | `bash code/03_evaluation/run_s2v8_predict.0.sh own_chromosome_configs/own_chromosome_predict.json` |
| `own_chromosome_configs/own_chromosome_predict.json` | Generated inference config matching all 231 final SSR fold models to their own held-out chromosomes. | Input to `deepcre_crosspredict.py`. | Adjust `EXTERNAL/` paths before re-running. |
| `crosspred_deployed/` (`job_crosspred_deployed.1.sh`, `deployed_models.tsv`, `configs/*.json`) | SGE runner + `26` per-(species × tissue) prediction JSONs scoring deployed SSR models across all 13 target genomes (13 leaf + 13 root configs; note `deployed_models.tsv` lists the 24 non-barley entries from the initial 12-species table). | Generates raw cross-species prediction CSVs behind Figure 2d–e, Supplementary Figure 1d–e, and cross-tissue evaluation. | Stage `EXTERNAL/environment` and `EXTERNAL/cluster_scratch`, then submit with `qsub`. |
| `crosspred_msr/` (`job_crosspred_msr.0.sh`, `{leaf,root}_models.txt`, `configs/*.json`) | SGE runner + `26` configs + exact `26` leaf / `26` root MSR model filename lists predicting all 52 MSR models across all 13 target genomes. | Generates raw MSR prediction CSVs behind Figure 2f–s and Supplementary Figure 1f–s. | Stage `EXTERNAL/environment` and `EXTERNAL/cluster_scratch`, then submit with `qsub`. |

---

## 2. Evaluation & Matrix Scripts (Python 3.11)

All scripts below read raw prediction CSVs (`EXTERNAL/predictions/`) and/or `data/`, producing the committed tables in [`results/evaluation/`](../../results/README.md):

| Script | Purpose | Manuscript outputs (`results/evaluation/`) | CLI Usage |
|---|---|---|---|
| `s2v8_validation_metrics.2.py` | Scores all 231 final SSR folds using 100 per-model CRC32-seeded balanced draws on `v8` paralog-filtered genes (`a`) and unbalanced metrics on all held-out chromosome genes (`b`), selecting the deployed fold per species × tissue by highest `auROC_v8`. | `final/S2v8_per_fold.tsv`, `final/S2v8_per_species.tsv`, `final/S2v8_per_tissue.tsv`, `final/deployed_models_{leaf,root}.tsv` | `python3 code/03_evaluation/s2v8_validation_metrics.2.py <final_models.tsv> <label_dir> <out_dir> <pred_dir> <v8.pickle> [<v8.pickle> ...]` |
| `supptable1_fold_counts.0.py` | Computes exact per-fold training vs. validation low/high/balanced gene counts alongside `v8` metrics. | `final/SuppTable1_per_fold_counts.tsv` (**Supplementary Table 1**) | `python3 code/03_evaluation/supptable1_fold_counts.0.py <S2v8_per_fold.tsv> <label_dir> <out.tsv> <pred_dir> [...]` |
| `crossspecies_matrices.2.py` | Builds 13 × 13 matrices of mean predicted probability per true class (`0` = low, `1` = high) across model species × target species. | `final/crossspecies_mean_{leaf,root}_class{0,1}.tsv` (**Figure 2d–e**, **Supp. Fig. 1d–e**) | `python3 code/03_evaluation/crossspecies_matrices.2.py <out_prefix> <label_dir> <tissue> <models.tsv> <pred_dir> [...]` |
| `msr_eval_dPm.2.py` | Evaluates MSR models (`all`, `dicots`, `monocots`) against held-out, seen, and outside-set targets (`100` balanced draws) vs. deployed SSR models. | `msr/msr_eval_{leaf,root}_{full,summary}.tsv` (**Figure 2f–s**, **Supp. Fig. 1f–s**) | `python3 code/03_evaluation/msr_eval_dPm.2.py <tissue> <msr_pred_dir> <ssr_pred_dir> <S2v8_per_species.tsv> <out_prefix>` |
| `msr_matrices.0.py` | Builds MSR 13 × 13 (and clade-specific 9 × 13 / 4 × 13) predicted-probability heatmap matrices. | `msr/suppfig_msr_mean_{leaf,root}_{all,dicots,monocots}_class{0,1}.tsv` | `python3 code/03_evaluation/msr_matrices.0.py <tissue> <msr_pred_dir> <out_prefix>` |
| `ssc_control_table.0.py` | Joins shuffled-sequence (`..._ctrlSSC.tsv`) and native/rescue training restart logs per fold. | `ssc_control/ssc_control_per_{fold,species}.tsv` (all 13 species incl. *H. vulgare*) | `python3 code/03_evaluation/ssc_control_table.0.py <restart_log_dir> <out_dir>` |
| `tissue_cross_eval.1.py` | Evaluates the 26 deployed SSR models across leaf vs. root expression labels restricted to genes scored in both tissues (shared across `s2v8` predictions). | `tissue_cross/tissue_cross_{full,own_species}.tsv` (**Supp. Fig. 2**) | `python3 code/03_evaluation/tissue_cross_eval.1.py <ssr_pred_dir> <out_prefix>` |
| `tissue_cross_msr.0.py` | Imports `tissue_cross_eval.1.py` and runs the matching leaf vs. root cross-tissue analysis for MSR models. | `tissue_cross/tissue_cross_msr_{full,switch_direction}.tsv` (**Supp. Fig. 2**) | `python3 code/03_evaluation/tissue_cross_msr.0.py <msr_pred_dir> <ssr_pred_dir> <out_prefix>` |
| `tissue_label_composition.1.py` | Tabulates genes into `stable`, `swapped`, `high to intermediate`, and `low to intermediate` across leaf & root. | `tissue_cross/tissue_label_composition.tsv` (**Supp. Fig. 2j**) | `python3 code/03_evaluation/tissue_label_composition.1.py <ssr_pred_dir> <out.tsv>` |
| `tissue_cross_switchers.1.py` | Scores leaf and root SSR models specifically on genes whose class flips between leaf and root. | `tissue_cross/tissue_switchers.tsv` | `python3 code/03_evaluation/tissue_cross_switchers.1.py <ssr_pred_dir> <out.tsv>` |
| `tissue_switch_direction.1.py` | Computes paired leaf vs. root probability contrasts (`P(leaf) - P(root)`) + 10,000-permutation tests on switched genes. | `tissue_cross/tissue_switch_direction.tsv`, `tissue_cross/tissue_switch_permutation.tsv` | `python3 code/03_evaluation/tissue_switch_direction.1.py <ssr_pred_dir> <out.tsv>` |
| `tissue_switch_performance.0.py` | Computes thresholded accuracy/MCC and threshold-free auROC on tissue-switching genes for SSR & MSR models. | `tissue_cross/tissue_switch_performance.tsv` (**Supp. Fig. 2m,o**) | `python3 code/03_evaluation/tissue_switch_performance.0.py <ssr_pred_dir> <msr_pred_dir> <out.tsv>` |
| `build_supplementary_data1.1.py` | Packages `final/`, `msr/`, and `ssc_control_per_fold.tsv` into the multi-sheet Excel workbook. *(Note: internal `LEGEND` string preserves working label `"Supplementary Data 6"` matching the original submission; published as Supplementary Data 1.)* | **Supplementary Data 1** `.xlsx` (offline-runnable directly from `results/evaluation`) | `python3 code/03_evaluation/build_supplementary_data1.1.py results/evaluation results/evaluation/ssc_control/ssc_control_per_fold.tsv <out.xlsx>` |
| `build_web_model_manifest.1.py` | Builds the 28-model web-tool deployment manifest (`26` SSR + `2` *B. oleracea* held-out `MSR_all` models) with SHA-256 verification. | `manifest.tsv` for `NAMlab/streamlit_deepcre_app` | `python3 code/03_evaluation/build_web_model_manifest.1.py <model_id_dir> <per_fold.tsv> <msr_dir> <prev_manifest|NONE> <out.tsv>` |

---

## How to Change or Extend

* **Add a new species to evaluation matrices and tissue analyses:** Update `ORDER` / `SP` / `LABELS` dictionaries at the top of `s2v8_validation_metrics.2.py`, `crossspecies_matrices.2.py`, `msr_eval_dPm.2.py`, `msr_matrices.0.py`, and `tissue_cross_eval.1.py`.
* **Change deployed MSR hold-out species on the web tool:** Update `DEPLOY_HELD_OUT` in `build_web_model_manifest.1.py` and `build_supplementary_data1.1.py`.
* **Smoke-test Supplementary Data 1 generation offline:** Run the `build_supplementary_data1.1.py` command in the table above using `results/evaluation`.
