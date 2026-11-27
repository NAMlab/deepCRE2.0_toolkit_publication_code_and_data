#!/bin/bash
REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"   # repository root
# Purpose: SGE CPU job evaluating the 26 deployed SSR models across target species (`deepcre_crosspredict.py`).
# Manuscript (Peleke et al., 2026): Generates cross-species prediction CSVs (`EXTERNAL/predictions/`) for Figure 2d–e & Supp. Fig. 1d–e.
# Usage: Stage EXTERNAL/environment and EXTERNAL/cluster_scratch, then `qsub code/03_evaluation/crosspred_deployed/job_crosspred_deployed.1.sh`.
# How to change or extend: Cap TF/OpenMP threads as shown to avoid CPU thread-allocation errors and update config paths.
#$ -N crosspred_dPm_r
#$ -tc 6
#$ -o $REPO/EXTERNAL/cluster_scratch/logs
#$ -wd $REPO/code/06_benchmarking/workflows/deepcre/src/deepCRE
#$ -j y
#$ -q Default.q
#$ -pe smp 8
#$ -l h_vmem=16g
source $REPO/EXTERNAL/environment/conda.sh
conda activate $REPO/EXTERNAL/environment/deepCRE_p38tf210_gpu
export PYTHONPATH=$REPO/code/06_benchmarking/workflows/deepcre/src CUDA_VISIBLE_DEVICES=""
export TF_NUM_INTRAOP_THREADS=4 TF_NUM_INTEROP_THREADS=2 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
CFG=$(ls "$REPO"/code/03_evaluation/crosspred_deployed/configs/*.json | sort | sed -n "${SGE_TASK_ID}p")
echo "task ${SGE_TASK_ID}: ${CFG}"
python deepcre_crosspredict.py -i "${CFG}"
