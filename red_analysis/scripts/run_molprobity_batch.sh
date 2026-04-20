#!/bin/bash
# Run MolProbity analysis for a list of targets sequentially
# Usage: bash run_molprobity_batch.sh <pipeline> <target1> <target2> ...

PYTHON=/home/agarwm5/miniconda3/envs/red_analysis/bin/python
SCRIPT=/data/p_csb_meiler/agarwm5/red_analysis/scripts/compute_rosetta_molprobity.py
PIPELINE=$1
shift

for TARGET in "$@"; do
    echo "=========================================="
    echo "Starting: $TARGET ($(date))"
    echo "=========================================="
    $PYTHON $SCRIPT $PIPELINE $TARGET 2>&1
    EXIT_CODE=$?
    echo "Finished: $TARGET (exit code: $EXIT_CODE, $(date))"
    echo ""
done

echo "=========================================="
echo "ALL DONE ($(date))"
echo "=========================================="
