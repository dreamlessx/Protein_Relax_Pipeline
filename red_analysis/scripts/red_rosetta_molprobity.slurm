#!/bin/bash
#SBATCH --job-name=red_ros_mp
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --array=1-257%30
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/red_analysis/logs/rosetta_mp_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/red_analysis/logs/rosetta_mp_%a.err

###############################################################################
# red_rosetta_molprobity.slurm - Phase 4: MolProbity for Rosetta structures
#
# COMMENT: This is the biggest compute job in Red's pipeline.
# Each target can have up to 780 .pdb.gz files × 2 pipelines = 1560 structures.
# At ~15-30 sec per structure, worst case is ~6-12 hours per target.
# 12-hour limit should be sufficient, but some targets may need re-runs.
#
# Only submit this AFTER Rosetta jobs are near-complete.
# Handles missing files gracefully - partial runs are fine.
#
# Concurrency reduced to 30 (from 50) to be a better cluster neighbor.
###############################################################################

set -euo pipefail

eval "$(~/miniconda3/bin/conda shell.bash hook 2>/dev/null)"
conda activate red_analysis

export MMTBX_CCP4_MONOMER_LIB="$CONDA_PREFIX/share/ccp4_mon_lib"

SCRIPT_DIR="/data/p_csb_meiler/agarwm5/red_analysis/scripts"

echo "=== Red Rosetta MolProbity: Task ${SLURM_ARRAY_TASK_ID} ==="
echo "Start: $(date)"

python3 "${SCRIPT_DIR}/compute_rosetta_molprobity.py" blue
python3 "${SCRIPT_DIR}/compute_rosetta_molprobity.py" green

echo "End: $(date)"
