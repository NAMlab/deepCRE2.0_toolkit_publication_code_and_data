#!/bin/bash
REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"   # repository root
# Purpose: SGE GPU array job (`1-6`) training leave-one-species-out MSR models (`all`, `dicots`, `monocots` x leaf/root).
# Manuscript (Peleke et al., 2026): Primary training run for all 52 multi-species models (`code/02_training/msr/configs/`).
# Usage: Stage EXTERNAL/environment and EXTERNAL/cluster_scratch, then `qsub code/02_training/msr/job_retrain_msr_dPmean.0.sh`.
# How to change or extend: Adjust `CFG` glob, conda activation, or `--seed` / `--max_attempts` flags.
#$ -N retrain_msr_dPm
#$ -t 1-6
#$ -o $REPO/EXTERNAL/cluster_scratch/logs
#$ -wd $REPO/code/06_benchmarking/workflows/deepcre/src/deepCRE
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=256g,gpu=1
#$ -P GPU_express.p
source $REPO/EXTERNAL/environment/conda.sh
conda activate $REPO/EXTERNAL/environment/deepCRE_p38tf210_gpu
CFG=$(ls "$REPO"/code/02_training/msr/configs/*.json | sort | sed -n "${SGE_TASK_ID}p")
echo "task ${SGE_TASK_ID}: ${CFG}"
python train_models_seed_restart.1.py --input "${CFG}" --seed 2026 --max_attempts 5 --plateau_loss 0.68
