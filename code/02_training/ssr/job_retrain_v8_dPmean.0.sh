#!/bin/bash
REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"   # repository root
# Purpose: SGE GPU array job (`1-26`) training single-species chromosome folds for all 13 species x leaf/root.
# Manuscript (Peleke et al., 2026): Primary SSR training run behind Figure 2, Supplementary Figure 1, and web-tool models.
# Usage: Stage EXTERNAL/environment and EXTERNAL/cluster_scratch, then `qsub code/02_training/ssr/job_retrain_v8_dPmean.0.sh`.
# How to change or extend: Replace conda env or JSON config directory (`CFG`) and adjust SGE `#$` directives for your cluster.
#$ -N retrain_v8_dPm
#$ -t 1-26
#$ -tc 6
#$ -o $REPO/EXTERNAL/cluster_scratch/logs
#$ -wd $REPO/deepCRE_reimplemented/src/deepCRE
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=160g,gpu=1
#$ -P GPU_express.p
source $REPO/EXTERNAL/environment/conda.sh
conda activate $REPO/EXTERNAL/environment/deepCRE_p38tf210_gpu
CFG=$(ls "$REPO"/code/02_training/ssr/configs/*.json | sort | sed -n "${SGE_TASK_ID}p")
echo "task ${SGE_TASK_ID}: ${CFG}"
python train_models_seed_restart.1.py --input "${CFG}" --seed 2026 --max_attempts 5 --plateau_loss 0.68
