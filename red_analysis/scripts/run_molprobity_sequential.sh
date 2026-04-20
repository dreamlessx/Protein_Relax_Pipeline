#!/bin/bash
# Run MolProbity analysis for Blue targets sequentially
# Each target gets its own Python process to avoid memory buildup
# Runs on login node since SLURM queue is congested

PYTHON=/home/agarwm5/miniconda3/envs/red_analysis/bin/python
SCRIPT=/data/p_csb_meiler/agarwm5/red_analysis/scripts/compute_rosetta_molprobity.py
LOGFILE=/data/p_csb_meiler/agarwm5/red_analysis/molprobity_sequential.log

TARGETS="1JTD 1JTG 1JWH 1JZD 1K5D 1K74 1KAC 1KKL 1KLU 1KTZ 1KXP 1LFD"

echo "Starting sequential MolProbity run ($(date))" > "$LOGFILE"

for TARGET in $TARGETS; do
    echo "=========================================="  | tee -a "$LOGFILE"
    echo "Starting: $TARGET ($(date))"                 | tee -a "$LOGFILE"
    echo "=========================================="  | tee -a "$LOGFILE"

    # Run in a subshell with ulimit to prevent runaway memory
    $PYTHON $SCRIPT blue $TARGET >> "$LOGFILE" 2>&1
    EXIT_CODE=$?

    echo "Finished: $TARGET (exit code: $EXIT_CODE, $(date))" | tee -a "$LOGFILE"

    # Check output file
    OUTFILE="/data/p_csb_meiler/agarwm5/red_analysis/metrics/rosetta_molprobity_blue_${TARGET}.tsv"
    if [ -f "$OUTFILE" ]; then
        ROWS=$(wc -l < "$OUTFILE")
        echo "  Output: $ROWS lines" | tee -a "$LOGFILE"
    else
        echo "  Output: NOT FOUND" | tee -a "$LOGFILE"
    fi
    echo "" | tee -a "$LOGFILE"
done

echo "=========================================="  | tee -a "$LOGFILE"
echo "ALL DONE ($(date))"                          | tee -a "$LOGFILE"
echo "=========================================="  | tee -a "$LOGFILE"
