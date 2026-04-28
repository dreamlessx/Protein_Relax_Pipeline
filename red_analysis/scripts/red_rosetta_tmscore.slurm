#!/bin/bash
#SBATCH --job-name=red_ros_tm
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --array=1-257%50
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/red_analysis/logs/rosetta_tm_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/red_analysis/logs/rosetta_tm_%a.err

###############################################################################
# red_rosetta_tmscore.slurm - TM-score for all Rosetta-relaxed structures
#
# COMMENT: Only submit this when Rosetta jobs are complete (or near-complete).
# Each target may have up to 780 pdb.gz files to score. With TMscore taking
# ~0.5s per file, worst case is ~7 minutes per target. 4 hour limit is generous.
#
# Run for partial data too - the script handles missing files gracefully.
###############################################################################

set -euo pipefail

SCRIPT_DIR="/data/p_csb_meiler/agarwm5/red_analysis/scripts"

echo "=== Red Rosetta TM-score: Task ${SLURM_ARRAY_TASK_ID} ==="
echo "Start: $(date)"

bash "${SCRIPT_DIR}/compute_rosetta_tmscore.sh" blue
bash "${SCRIPT_DIR}/compute_rosetta_tmscore.sh" green

echo "End: $(date)"
