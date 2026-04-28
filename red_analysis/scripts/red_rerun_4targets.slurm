#!/bin/bash
#SBATCH --job-name=red_mp_rerun
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --array=1-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/red_analysis/logs/rerun_mp_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/red_analysis/logs/rerun_mp_%a.err

###############################################################################
# red_rerun_4targets.slurm — Re-run Phase 4 MolProbity for 3 Blue targets
# (1A2K already done, so only 1ATN, 1BKD, 1D6R need processing)
#
# Array task mapping:
#   1 = 1ATN
#   2 = 1BKD
#   3 = 1D6R
#
# Pre-run cleanup: deletes any existing output file with NAs or header-only
###############################################################################

set -euo pipefail

eval "$(~/miniconda3/bin/conda shell.bash hook 2>/dev/null)"
conda activate red_analysis

export MMTBX_CCP4_MONOMER_LIB="$CONDA_PREFIX/share/ccp4_mon_lib"

SCRIPT_DIR="/data/p_csb_meiler/agarwm5/red_analysis/scripts"
METRICS_DIR="/data/p_csb_meiler/agarwm5/red_analysis/metrics"

# Map array task ID to target
declare -a TARGETS=("1ATN" "1BKD" "1D6R")
TARGET=${TARGETS[$((SLURM_ARRAY_TASK_ID - 1))]}

echo "=== Red MolProbity Rerun: ${TARGET} (Task ${SLURM_ARRAY_TASK_ID}) ==="
echo "Start: $(date)"

OUTFILE="${METRICS_DIR}/rosetta_molprobity_blue_${TARGET}.tsv"

# Check if file exists and needs cleanup
if [ -f "$OUTFILE" ]; then
    NA_COUNT=$(awk -F'\t' 'NR>1 && $6=="NA" {count++} END {print count+0}' "$OUTFILE")
    TOTAL_ROWS=$(awk -F'\t' 'NR>1 {count++} END {print count+0}' "$OUTFILE")
    echo "Existing file: $TOTAL_ROWS rows, $NA_COUNT NA clashscores"

    if [ "$NA_COUNT" -gt 0 ] || [ "$TOTAL_ROWS" -eq 0 ]; then
        echo "Deleting file (has NAs or is empty) and re-running..."
        rm -f "$OUTFILE"
    else
        echo "SKIP: File is clean with $TOTAL_ROWS rows"
        exit 0
    fi
fi

python3 "${SCRIPT_DIR}/compute_rosetta_molprobity.py" blue "${TARGET}"

# Verify output
if [ -f "$OUTFILE" ]; then
    ROWS=$(awk -F'\t' 'NR>1 {count++} END {print count+0}' "$OUTFILE")
    NAS=$(awk -F'\t' 'NR>1 && $6=="NA" {count++} END {print count+0}' "$OUTFILE")
    echo "Result: $ROWS rows, $NAS NA clashscores"
else
    echo "WARNING: Output file not created!"
fi

echo "End: $(date)"
