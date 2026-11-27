#!/bin/bash
REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"   # repository root
# Purpose: SGE GPU job resuming interrupted all-species MSR tasks with thread caps via `--only`.
# Manuscript (Peleke et al., 2026): Completed held-out leaf (`Zmay`, `Hvul`) and root (`Dcar`–`Hvul`)
#   multi-species models.
# Usage: Stage EXTERNAL/environment and EXTERNAL/cluster_scratch, then `qsub code/02_training/msr/job_finish_msr_all_dPmean.0.sh` (or run `train_models_seed_restart.2.py` locally).
# How to change or extend: Edit `ONLY`, `CFG`, conda env, and SGE `#$` resource directives for your cluster.
#$ -N finish_msr_dPm
#$ -t 1-2
#$ -o $REPO/EXTERNAL/cluster_scratch/logs
#$ -wd $REPO/code/06_benchmarking/workflows/deepcre/src/deepCRE
#$ -j y
#$ -pe smp 1
# Note: h_vmem is requested per SGE slot
#$ -l h_vmem=300g,gpu=1
#$ -P GPU_express.p
source $REPO/EXTERNAL/environment/conda.sh
conda activate $REPO/EXTERNAL/environment/deepCRE_p38tf210_gpu
export TF_NUM_INTRAOP_THREADS=4 TF_NUM_INTEROP_THREADS=2 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
if [ "${SGE_TASK_ID}" -eq 1 ]; then
  CFG="$REPO/code/02_training/msr/configs/01_MSR_all_leaf.json"; ONLY=Zmay_dPm_leaf,Hvul_dPm_leaf
else
  CFG="$REPO/code/02_training/msr/configs/02_MSR_all_root.json"; ONLY=Dcar_dPm_root,Gmax_dPm_root,Slyc_dPm_root,Osat_dPm_root,Sbic_dPm_root,Zmay_dPm_root,Hvul_dPm_root
fi
echo "task ${SGE_TASK_ID}: ${CFG} only ${ONLY}"
python train_models_seed_restart.2.py --input "${CFG}" --seed 2026 --max_attempts 5 --plateau_loss 0.68 --only "${ONLY}"
