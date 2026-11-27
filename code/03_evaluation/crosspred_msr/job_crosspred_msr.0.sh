#!/bin/bash
REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"   # repository root
# Purpose: SGE CPU array job (`1-26`) running `deepcre_crosspredict.py` for all 52 MSR models across 13 target genomes.
# Manuscript (Peleke et al., 2026): Generates MSR cross-prediction CSVs (`EXTERNAL/predictions/`) for Figure 2f–s & Supp. Fig. 1f–s.
# Usage: Stage EXTERNAL/environment and EXTERNAL/cluster_scratch, then `qsub code/03_evaluation/crosspred_msr/job_crosspred_msr.0.sh`.
# How to change or extend: Update config folder paths and `-t` array size when adding new MSR species.
#$ -t 1-26
#$ -N crosspred_msr
#$ -o $REPO/EXTERNAL/cluster_scratch/logs
#$ -wd $REPO/code/06_benchmarking/workflows/deepcre/src/deepCRE
#$ -j y
#$ -q Default.q
#$ -pe smp 4
#$ -l h_vmem=24g
source $REPO/EXTERNAL/environment/conda.sh
conda activate $REPO/EXTERNAL/environment/deepCRE_p38tf210_gpu
export PYTHONPATH=$REPO/code/06_benchmarking/workflows/deepcre/src CUDA_VISIBLE_DEVICES=""
export TF_NUM_INTRAOP_THREADS=4 TF_NUM_INTEROP_THREADS=2 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
CFG=$(ls "$REPO"/code/03_evaluation/crosspred_msr/configs/*.json | sort | sed -n "${SGE_TASK_ID}p")
echo "task ${SGE_TASK_ID}: ${CFG}"
python deepcre_crosspredict.py -i "${CFG}"
