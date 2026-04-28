#!/bin/bash --norc
#SBATCH --account=csb_gpu_acc
#SBATCH --constraint=csbtmp
#SBATCH --partition=batch_gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:nvidia_rtx_a6000:1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --job-name=green_axtal_fill
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/logs/green_axtal_fill_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/logs/green_axtal_fill_%A_%a.err
#SBATCH --array=1,2,6,7,16,20,21,22,24,25,26,27,28,30,40,42,44,45,46,47,48,53,56,57,67,69,70,76,83,86,89,97,102,105,106,107,108,111,115,116,117,121,126,130,131,133,135,136,137,140,141,145,149,156,158,159,161,162,169,170,173,174,175,176,177,178,179,180,181,183,185,186,187,188,189,190,192,198,200,203,206,209,215,216,218,225,226,228,229,230,237,240,249,251,254
#SBATCH --exclude=cn1340

# AMBER crystal relaxation for 95 Green targets missing it
# Uses af_dirlist.txt indexing to match Green's data layout

set -eo pipefail

# ACCRE migrated module system - load new stack first
setup_accre_software_stack 2>/dev/null || true
source /sb/apps/alphafold232/miniconda3/bin/activate af232

set -u

BASE="/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking"
DIRLIST="${BASE}/af_dirlist.txt"
CLEANED="${BASE}/cleaned"
SCRIPT="/data/p_csb_meiler/agarwm5/protein_pipeline/scripts/amber_relax_crystal_v3.py"

D=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$DIRLIST")
[[ -z "$D" || ! -d "$D" ]] && { echo "No target for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

TARGET=$(basename "$D")
CRYSTAL="$CLEANED/${TARGET}.pdb"

echo "AMBER crystal fill (Green): $TARGET on $(hostname)"
nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || true

python "$SCRIPT" "$D" "$CRYSTAL"
