# Trained deepCRE Models (`models/`)

All trained Keras `.h5` weights (~4.8 MB each) are stored in this repository:

| Directory | Models |
|---|---|
| `ssr/` | 232 single-species models (SSR) |
| `ssc/` | 217 shuffled-sequence controls (SSC) |
| `msr/` | 52 multi-species leave-one-species-out models (MSR) |
| `peleke2024/` | 45 models from Peleke et al. (2024) |

Files carry the standardized names (`new_id`) of `model_id_mapping.tsv`; the configs in `code/` point to `models/<group>/<new_id>`.

---

## 2026 Model Inventory (`model_id_mapping.tsv`)

`model_id_mapping.tsv` indexes **501** archived 2026 `.h5` files (trained with `deepCRE_reimplemented` at commit `19f6db7` via `code/02_training/`):

| Former release archive | Model group | Count in `model_id_mapping.tsv` | Notes & Manuscript Role |
|---|---|---|---|
| `deepCRE_TK_2026_models_SSR.tar` | Single-species (`single-species`) | **232** (`231` `in_final_set=True`: `115` leaf + `116` root; `1` `superseded/partial` *C. sativa* leaf fold) | One CNN per held-out chromosome across 13 species × leaf/root. The highest-validation-auROC converged fold per species × tissue (**26** models, marked `deployed=True`) is deployed on the web tool. |
| `deepCRE_TK_2026_models_SSC.tar` | Shuffled-sequence controls (`SSC control`) | **217** (`12` species × leaf/root) | Identical chromosome folds and hyperparameters with promoter/terminator nucleotides randomly shuffled during training (`model_case="ssc"`). Note: *H. vulgare* SSC `.h5` weights are not packaged in this tar archive, although *H. vulgare* SSC training metrics (14 folds) are included in `results/evaluation/ssc_control/`. |
| `deepCRE_TK_2026_models_MSR.tar` | Multi-species leave-one-species-out (`multi-species`) | **52** (`26` all-13-species, `18` 9-dicots-only, `8` 4-monocots-only across leaf & root) | One model per held-out species × tissue. Two all-species models held out on *B. oleracea* (`MSR_all_M0X0.75dPm_Bole_...` leaf and `MSR_all_M1X0.75dPm_Bole_...` root) are selected for web-tool deployment by `code/03_evaluation/build_web_model_manifest.1.py` (making **28** deployed 2026 models in total). |

### `model_id_mapping.tsv` Columns

| Column | Meaning |
|---|---|
| `kind` | `single-species`, `SSC control`, `multi-species (all)`, `multi-species (dicots)`, or `multi-species (monocots)` |
| `code`, `tissue` | 4-letter species code (`Atha`…`Zmay`) and tissue (`leaf`, `root`); for MSR, the held-out species |
| `fold_or_heldout` | Held-out chromosome accession (`NC_*`/`NW_*` for SSR/SSC) or held-out species code (MSR) |
| `old_file` | Original filename written during training (`code/02_training/`) and referenced in `results/evaluation/` tables |
| `new_id` | Standardized publication filename stored in `models/ssr/`, `models/ssc/` or `models/msr/` |
| `in_final_set` | `True` for the 231 converged SSR folds and 52 MSR models evaluated in the manuscript |
| `status` | `converged`, `converged_rescue`, `control`, or `superseded/partial` |
| `deployed` | `True` for the 26 deployed SSR fold models |

---

## Peleke et al. (2024) Models (`peleke2024_models.tsv`)

`models/peleke2024/` holds **45** models with SHA-256 hashes indexed in `peleke2024_models.tsv`:
* `per_chromosome_2024` (**37**): per-chromosome hold-out models for *A. thaliana* (5), *S. lycopersicum* (12), *S. bicolor* (10), and *Z. mays* (10).
* `web_tool_2024` (**8**): leaf and root models of those four species served on the web tool as "Peleke et al., 2024" and compared in `code/05_supplementary_figures/` (Supplementary Figures 5 and 6).

---

## Swapping or Extending Models

To add or replace models in downstream scripts:
1. Place new `.h5` weights under `models/<group>/` and regenerate `model_id_mapping.tsv` / `final_models.tsv` using `code/02_training/s2v8_select_models.1.py` and `code/02_training/assign_model_ids.0.py`.
2. Regenerate the web-tool manifest with `code/03_evaluation/build_web_model_manifest.1.py` (update `DEPLOY_HELD_OUT` in that script to change which MSR hold-out species is served).
