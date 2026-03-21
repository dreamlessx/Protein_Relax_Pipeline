#!/bin/bash
#SBATCH --job-name=red_mp4_12
#SBATCH --partition=batch
#SBATCH --account=p_csb_meiler
#SBATCH --array=1-12
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/red_analysis/logs/mp4_%a_%A.log
#SBATCH --error=/data/p_csb_meiler/agarwm5/red_analysis/logs/mp4_%a_%A.err

# Phase 4 MolProbity for 12 Blue targets
# Array index maps to target

TARGETS=(1GHQ 1GL1 1GLA 1GP2 1GPW 1GRN 1GXD 1H1V 1H9D 1HCF 1HE1 1HE8)
TARGET=${TARGETS[$SLURM_ARRAY_TASK_ID-1]}

echo "=== Phase 4 MolProbity: ${TARGET} (blue) ==="
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Started: $(date)"

/home/agarwm5/miniconda3/envs/red_analysis/bin/python -u \
    /data/p_csb_meiler/agarwm5/red_analysis/scripts/compute_rosetta_molprobity.py \
    blue ${TARGET}

echo "Finished: $(date)"
