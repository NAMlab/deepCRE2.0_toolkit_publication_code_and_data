#!/bin/bash
REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"   # repository root
# Purpose: SGE GPU array job (`1-26`) training single-attempt shuffled-sequence control (`ssc`) models.
# Manuscript (Peleke et al., 2026): Shuffled-input negative controls across 13 species x leaf/root (`ssc_control/configs/`).
# Usage: Stage EXTERNAL/environment and EXTERNAL/cluster_scratch, then `qsub code/02_training/ssc_control/job_ssc_control_dPmean.0.sh`.
# How to change or extend: Point `CFG` to your JSON config folder or adjust thread caps / memory limits.
#$ -N ssc_ctrl_dPm
#$ -t 1-26
#$ -tc 6
#$ -o $REPO/EXTERNAL/cluster_scratch/logs
#$ -wd $REPO/code/06_benchmarking/workflows/deepcre/src/deepCRE
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=160g,gpu=1
#$ -P GPU_express.p
source $REPO/EXTERNAL/environment/conda.sh
conda activate $REPO/EXTERNAL/environment/deepCRE_p38tf210_gpu
export TF_NUM_INTRAOP_THREADS=4 TF_NUM_INTEROP_THREADS=2 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
CFG=$(ls "$REPO"/code/02_training/ssc_control/configs/*.json | sort | sed -n "${SGE_TASK_ID}p")
echo "task ${SGE_TASK_ID}: ${CFG}"
python train_models_seed_restart.2.py --input "${CFG}" --seed 2026 --max_attempts 1 --plateau_loss 0.68
