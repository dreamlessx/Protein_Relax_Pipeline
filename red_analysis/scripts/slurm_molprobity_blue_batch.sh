#!/bin/bash
#SBATCH --job-name=red_mp_blue
#SBATCH --partition=batch
#SBATCH --account=p_csb_meiler
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/red_analysis/logs/molprobity_blue_batch_%j.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/red_analysis/logs/molprobity_blue_batch_%j.err

# Red Analysis Pipeline — Phase 4 MolProbity for 12 Blue targets
# Run sequentially on a compute node to avoid login node OOM

PYTHON=/home/agarwm5/miniconda3/envs/red_analysis/bin/python
SCRIPT=/data/p_csb_meiler/agarwm5/red_analysis/scripts/compute_rosetta_molprobity.py

TARGETS="1JTD 1JTG 1JWH 1JZD 1K5D 1K74 1KAC 1KKL 1KLU 1KTZ 1KXP 1LFD"

for TARGET in $TARGETS; do
    echo "=========================================="
    echo "Starting: $TARGET ($(date))"
    echo "=========================================="
    $PYTHON $SCRIPT blue $TARGET 2>&1
    EXIT_CODE=$?
    echo "Finished: $TARGET (exit code: $EXIT_CODE, $(date))"
    echo ""
done

echo "=========================================="
echo "ALL DONE ($(date))"
echo "=========================================="

# Summary
echo ""
echo "=== SUMMARY ==="
for TARGET in $TARGETS; do
    FILE=/data/p_csb_meiler/agarwm5/red_analysis/metrics/rosetta_molprobity_blue_${TARGET}.tsv
    if [ -f "$FILE" ]; then
        ROWS=$(wc -l < "$FILE")
        echo "$TARGET: $ROWS lines ($(($ROWS - 1)) data rows)"
    else
        echo "$TARGET: NOT FOUND"
    fi
done
