#!/bin/bash
#SBATCH --job-name=red_tmscore
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --array=1-257%100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/red_analysis/logs/tmscore_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/red_analysis/logs/tmscore_%a.err

###############################################################################
# red_tmscore.slurm - Compute TM-score, RMSD, GDT for all pre-Rosetta
# structures against crystal references.
#
# Runs BOTH Blue and Green pipelines per target (one SLURM task per target).
# Pre-Rosetta only: AF, Boltz, AMBER (all 257/257 complete for both pipelines).
#
# Each target produces ~50 metric rows:
#   5 sources × 5 models × 2 pipelines = 50 TMscore computations
#   (each takes <1 sec, so 1 hour is very generous)
###############################################################################

set -euo pipefail

SCRIPT_DIR="/data/p_csb_meiler/agarwm5/red_analysis/scripts"

echo "=== Red TM-score: Task ${SLURM_ARRAY_TASK_ID} ==="
echo "Start: $(date)"

# Run for Blue pipeline
bash "${SCRIPT_DIR}/compute_tmscore.sh" blue

# Run for Green pipeline
bash "${SCRIPT_DIR}/compute_tmscore.sh" green

echo "End: $(date)"
