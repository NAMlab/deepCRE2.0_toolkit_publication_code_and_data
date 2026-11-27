# Model Training (`code/02_training/`)

Wrappers, SGE cluster job scripts, JSON configurations, and model selection/naming utilities used to train all 2026 deepCRE2.0 models reported in Peleke et al. (2026).

---

## Core Architecture & Wrapper (`train_models_seed_restart.{1,2}.py`)

All models were trained using the unmodified CNN architecture and data pipeline in **`code/06_benchmarking/workflows/deepcre/src/deepCRE/train_models.py`** (commit `19f6db7` of https://github.com/NAMlab/deepCRE2.0, formerly `deepCRE_reimplemented`; the submodule is pinned at `a37db45`, which differs only by skipping genes whose flanks reach past a contig end instead of stopping with an error; 3,020 bp one-hot gene windows `[1,000 bp gUR + 500 bp gTUR + 20 bp N-pad + 500 bp gTDR + 1,000 bp gDR]`, 3 Conv1D+MaxPool+Dropout blocks [64/128/64 filters, kernel 8, pool 8, dropout 0.25], Dense `128 → 64 → 1` sigmoid, Adam `lr = 1e-4`, binary cross-entropy, batch `64`, up to `100` epochs with EarlyStopping patience `10` and ReduceLROnPlateau patience `5`).

Because `train_models.py` did not fix random seeds and occasional chromosome folds stalled near constant `0.5` output (`val_loss ≈ ln(2) = 0.693`), the lightweight runtime wrapper **`train_models_seed_restart.{1,2}.py`** (placed next to `train_models.py` during execution) patches two functions in memory without editing `train_models.py`:

* `train_deep_cre`: sets `random`, `numpy`, and `tensorflow` seeds (`--seed 2026`) before each fold so `balance_dataset()`, train/validation splits, and SSC nucleotide shuffling are deterministic.
* `deep_cre`: runs up to `--max_attempts` weight re-initializations (`seed + attempt`, `attempt = first_attempt .. max_attempts - 1`) until best-checkpoint `val_loss <= --plateau_loss` (`0.68`), logging every attempt to `restart_log_<output_name>.tsv`.
* `.1.py` adds `--first_attempt` (used for rescue seeds `2031`–`2045`); `.2.py` adds `--only A,B` to run specific interrupted entities (MSR/SSC/TAIR10).

---

## Training Regimes, Cluster Scripts, and Configs

All JSON config files under `ssr/configs/`, `ssr/configs_rescue/`, `msr/configs/`, `ssc_control/configs/`, and `tair10_exception/` use paths relative to the repository root (`data/...`). SGE `.sh` job scripts resolve `$REPO` dynamically and reference `$REPO/EXTERNAL/environment/` and `$REPO/EXTERNAL/cluster_scratch/` (see `EXTERNAL_INPUTS.md`).

| Subfolder | `model_case` & Purpose | Cluster script | Config count |
|---|---|---|---|
| `ssr/` | Single-species per-chromosome cross-validation across 13 species × leaf/root (`configs/`), plus rescue run (`configs_rescue/` = *C. sativa* leaf) for folds stalling at `seed 2026`–`2030`. | `ssr/job_retrain_v8_dPmean.0.sh` (`26` array tasks, `--seed 2026 --max_attempts 5`); `ssr/job_rescue_v8_dPmean.0.sh` (`--first_attempt 5 --max_attempts 15`) | `26` base + `1` rescue JSON |
| `msr/` | Multi-species leave-one-species-out models across all 13 species (`01`/`02`), 9 dicots (`03`/`04`), and 4 monocots (`05`/`06`) in leaf and root. | `msr/job_retrain_msr_dPmean.0.sh` (`6` array tasks); `msr/job_finish_msr_all_dPmean.0.sh` (finishes held-out species via `--only`) | `6` JSONs |
| `ssc_control/` | Shuffled-sequence control (`ssc`): promoter and terminator nucleotides of training genes randomly shuffled (`np.random.shuffle`), validation sequences unshuffled, single attempt (`--max_attempts 1`). | `ssc_control/job_ssc_control_dPmean.0.sh` (`26` array tasks, `13` species × leaf/root) | `26` JSONs |
| `tair10_exception/` | *A. thaliana* models trained on curated TAIR10 RefSeq gene windows instead of Helixer windows. See [`tair10_exception/README.md`](tair10_exception/README.md). | `tair10_exception/job_train_atha_tair10.0.sh` (`2` array tasks, `--seed 2026 --max_attempts 5`) | `2` JSONs |

---

## Post-Training Selection and Model Naming Scripts

| File | Purpose | Manuscript role | Usage |
|---|---|---|---|
| `s2v8_select_models.1.py` | Reads trained `.h5` filenames and `restart_log_*.tsv` files against expected chromosome folds (`S2_fold_sizes` / `SuppTable1` list) and resolves final converged vs. rescue vs. superseded models. | Generates `final_models.tsv` (`231` converged single-species fold entries across `13` species × `2` tissues). | `python3 code/02_training/s2v8_select_models.1.py <model_dir> <restart_log_dir> <expected_folds.tsv> <out.tsv>` |
| `assign_model_ids.0.py` | Assigns standardized publication model IDs (`S0X`/`S1X`/`M0X`/`M1X`/`M0M`/`M1M` + `0.75dPm`), creates hard links, and emits the mapping index. | Generates [`models/model_id_mapping.tsv`](../../models/model_id_mapping.tsv) (`501` archived 2026 `.h5` entries). | `python3 code/02_training/assign_model_ids.0.py <ssr_model_dir> <msr_model_dir> <final_models.tsv> <S2v8_per_fold.tsv> <out_dir>` |
| `final_models.tsv` | Resolved table of the 231 final SSR fold models (`tissue`, `code`, `chrom`, `file`, `status`, `seed`, validation metrics). | Input to `code/03_evaluation/s2v8_make_predict_config.1.py` and `s2v8_validation_metrics.2.py`. | — |

---

## How to Change or Extend

* **Train a new species or annotation:** Create a JSON modeled on `ssr/configs/01_Atha_leaf.json`, copy `train_models_seed_restart.2.py` into `code/06_benchmarking/workflows/deepcre/src/deepCRE/`, ensure `data/*/*.gz` are decompressed and `.fna` genomes exist, then run from `code/06_benchmarking/workflows/deepcre/src/deepCRE/`:
  ```bash
  python train_models_seed_restart.2.py --input path/to/config.json --seed 2026 --max_attempts 5 --plateau_loss 0.68
  ```
* **Use `--only` for partial restarts:** Pass `--only NC_025691.1` (for SSR fold chromosome) or `--only Hvul_dPm_leaf` (for MSR species name) to train only failed folds/species without overwriting completed ones.
