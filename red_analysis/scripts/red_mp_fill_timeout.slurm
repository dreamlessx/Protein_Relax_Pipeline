#!/bin/bash
#SBATCH --job-name=red_mp_to_fill
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/red_analysis/logs/red_mp_to_fill_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/red_analysis/logs/red_mp_to_fill_%a.err
#SBATCH --exclude=cn1340
#SBATCH --array=83,215
set -euo pipefail
eval "$(~/miniconda3/bin/conda shell.bash hook 2>/dev/null)"
conda activate red_analysis
export MMTBX_CCP4_MONOMER_LIB="$CONDA_PREFIX/share/ccp4_mon_lib"
SCRIPT_DIR="/data/p_csb_meiler/agarwm5/red_analysis/scripts"
echo "=== Red Rosetta MolProbity RETRY (48h): Task ${SLURM_ARRAY_TASK_ID} ==="
echo "Start: $(date)"
python3 "${SCRIPT_DIR}/compute_rosetta_molprobity.py" blue
python3 "${SCRIPT_DIR}/compute_rosetta_molprobity.py" green
echo "End: $(date)"
