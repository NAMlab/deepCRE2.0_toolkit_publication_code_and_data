#!/bin/bash
REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"   # repository root
# Purpose: Run `deepcre_crosspredict.py` locally with `EXTERNAL/environment/python` on a config JSON in `code/03_evaluation/`.
# Manuscript (Peleke et al., 2026): Generates own-chromosome `S2v8` validation predictions for `s2v8_validation_metrics.2.py`.
# Usage: bash code/03_evaluation/run_s2v8_predict.0.sh own_chromosome_configs/own_chromosome_predict.json
# How to change or extend: Point `EXTERNAL/environment/` to your Python 3.8 + TensorFlow 2.10 environment.
cd $REPO/deepCRE_reimplemented/src/deepCRE
export LD_LIBRARY_PATH=$REPO/EXTERNAL/environment/lib PYTHONPATH=$REPO/deepCRE_reimplemented/src/
W=$REPO/code/03_evaluation/
LOG=$W/${1%.json}.log
echo "$(date +%T) start $1"
$REPO/EXTERNAL/environment/python deepcre_crosspredict.py -i $W/$1 > $LOG 2>&1
echo "$(date +%T) end exit $?"
grep -E "passed\.|failed with error" $LOG | tail -3
