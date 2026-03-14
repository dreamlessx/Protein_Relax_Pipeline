#!/bin/bash
#SBATCH --job-name=red_molprob
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --array=1-257%50
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/red_analysis/logs/molprobity_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/red_analysis/logs/molprobity_%a.err

###############################################################################
# red_molprobity.slurm — MolProbity validation for all pre-Rosetta structures
#
# COMMENT: This is Phase 3 of Red analysis. Uses cctbx Python API with
# individual validators (ramalyze, rotalyze, cbetadev, clashscore) plus
# geometry restraints for RMS bonds/angles.
#
# Per target: 1 crystal + 5 af_relaxed + 5 af_unrelaxed + 5 boltz + 5 amber_af
#             + 5 amber_boltz = 26 structures × 2 pipelines = 52 MolProbity runs
#
# Each structure takes ~10-30 seconds. Worst case ~25 min/target.
# 4 hour limit is generous.
#
# REQUIRES: conda env red_analysis with cctbx-base, probe, reduce (wrapper)
#
# Key patches for this cluster:
#   1. reduce wrapper (forces return code 0 — reduce v4.16 returns 255 on success)
#   2. SS link alias (CCP4 lib uses 'disulf', cctbx expects 'SS')
#   3. probe/reduce module detection monkey-patch
#   4. H-stripping for AF-relaxed and AMBER structures (424 extra H atoms)
###############################################################################

set -euo pipefail

# Activate conda environment
eval "$(~/miniconda3/bin/conda shell.bash hook 2>/dev/null)"
conda activate red_analysis

# Set monomer library path
export MMTBX_CCP4_MONOMER_LIB="$CONDA_PREFIX/share/ccp4_mon_lib"

SCRIPT_DIR="/data/p_csb_meiler/agarwm5/red_analysis/scripts"

echo "=== Red MolProbity: Task ${SLURM_ARRAY_TASK_ID} ==="
echo "Start: $(date)"

python3 "${SCRIPT_DIR}/compute_molprobity_py.py" blue
python3 "${SCRIPT_DIR}/compute_molprobity_py.py" green

echo "End: $(date)"
