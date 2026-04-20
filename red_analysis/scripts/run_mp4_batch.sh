#!/bin/bash
# Run Phase 4 MolProbity for 12 Blue targets sequentially
# Each target is run one at a time

PYTHON=/home/agarwm5/miniconda3/envs/red_analysis/bin/python
SCRIPT=/data/p_csb_meiler/agarwm5/red_analysis/scripts/compute_rosetta_molprobity.py
LOGDIR=/data/p_csb_meiler/agarwm5/red_analysis/logs

TARGETS="1GHQ 1GL1 1GLA 1GP2 1GPW 1GRN 1GXD 1H1V 1H9D 1HCF 1HE1 1HE8"

echo "=== Phase 4 MolProbity Batch Run ==="
echo "Started: $(date)"
echo ""

for TARGET in $TARGETS; do
    echo ">>> Starting $TARGET at $(date)"
    $PYTHON -u $SCRIPT blue $TARGET 2>&1 | tee ${LOGDIR}/mp4_${TARGET}.log
    echo ">>> Finished $TARGET at $(date)"
    echo ""
done

echo "=== All 12 targets complete ==="
echo "Finished: $(date)"

# Summary
echo ""
echo "=== Output Summary ==="
for TARGET in $TARGETS; do
    OUTFILE=/data/p_csb_meiler/agarwm5/red_analysis/metrics/rosetta_molprobity_blue_${TARGET}.tsv
    if [ -f "$OUTFILE" ]; then
        LINES=$(wc -l < "$OUTFILE")
        DATAROWS=$((LINES - 1))
        NA_COUNT=$(grep -c "	NA	" "$OUTFILE" 2>/dev/null || echo 0)
        echo "$TARGET: $DATAROWS data rows, $NA_COUNT rows with NA"
    else
        echo "$TARGET: MISSING"
    fi
done

TOTAL=$(cat /data/p_csb_meiler/agarwm5/red_analysis/metrics/rosetta_molprobity_blue_1G*.tsv /data/p_csb_meiler/agarwm5/red_analysis/metrics/rosetta_molprobity_blue_1H*.tsv 2>/dev/null | grep -v "^target" | wc -l)
echo ""
echo "Total data rows across all 12 files: $TOTAL"
