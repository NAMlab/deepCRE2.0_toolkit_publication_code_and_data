#!/bin/bash
REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"   # repository root
# Purpose: SGE GPU array job re-initializing chromosome folds that remained on the 0.693 loss plateau after 5 seeds.
# Manuscript (Peleke et al., 2026): Rescue training using `train_models_seed_restart.1.py` (`--first_attempt 5 --max_attempts 15`).
# Usage: Stage EXTERNAL/environment and EXTERNAL/cluster_scratch, then `qsub code/02_training/ssr/job_rescue_v8_dPmean.0.sh`.
# How to change or extend: Place rescue JSONs in config folder and update `--first_attempt` / `--max_attempts`.
#$ -N rescue_v8_dPm
#$ -o $REPO/EXTERNAL/cluster_scratch/logs
#$ -wd $REPO/deepCRE_reimplemented/src/deepCRE
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=160g,gpu=1
#$ -P GPU_express.p
source $REPO/EXTERNAL/environment/conda.sh
conda activate $REPO/EXTERNAL/environment/deepCRE_p38tf210_gpu
CFG=$(ls "$REPO"/code/02_training/ssr/configs_rescue/*.json | sort | sed -n "${SGE_TASK_ID:-1}p")
echo "task ${SGE_TASK_ID}: ${CFG}"
python train_models_seed_restart.1.py --input "${CFG}" --seed 2026 --first_attempt 5 --max_attempts 15 --plateau_loss 0.68
