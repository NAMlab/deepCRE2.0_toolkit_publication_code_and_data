#!/bin/bash
REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"   # repository root
# Purpose: SGE GPU array job (`1-2`) training *A. thaliana* leaf & root models on curated TAIR10 RefSeq windows.
# Manuscript (Peleke et al., 2026): Trains companion `Atha_TAIR10_*` models for Supplementary Note 4 (Run F).
# Usage: Stage `EXTERNAL/environment` & `EXTERNAL/cluster_scratch`, then `qsub code/02_training/tair10_exception/job_train_atha_tair10.0.sh`.
# How to change or extend: Add configs under `code/02_training/tair10_exception/` or adjust SGE `#$` array range `-t`.
#$ -N atha_tair10
#$ -t 1-2
#$ -o $REPO/EXTERNAL/cluster_scratch/logs
#$ -wd $REPO/code/06_benchmarking/workflows/deepcre/src/deepCRE
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=64g,gpu=1
#$ -P GPU_express.p
source $REPO/EXTERNAL/environment/conda.sh
conda activate $REPO/EXTERNAL/environment/deepCRE_p38tf210_gpu
CFG=$(ls "$REPO"/code/02_training/tair10_exception/*.json | sort | sed -n "${SGE_TASK_ID}p")
echo "task ${SGE_TASK_ID}: ${CFG}"
python train_models_seed_restart.2.py --input "${CFG}" --seed 2026 --max_attempts 5 --plateau_loss 0.68
