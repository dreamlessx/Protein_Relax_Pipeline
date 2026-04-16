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
#SBATCH --job-name=blue_axtal_fill
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/blue_axtal_fill_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/blue_axtal_fill_%A_%a.err
#SBATCH --array=1,2,6,7,16,17,21,22,23,25,26,27,28,29,31,41,43,45,46,48,49,50,55,58,59,70,72,73,79,83,87,90,93,101,106,110,111,112,113,116,120,122,123,127,132,136,137,139,141,142,143,146,147,151,155,162,164,165,167,168,176,177,181,182,183,184,185,186,187,188,189,190,192,194,195,196,197,198,199,201,207,209,212,215,219,225,226,228,235,236,238,239,240,247,250
#SBATCH --exclude=cn1340

# AMBER crystal relaxation for 95 targets missing it
# GPU job — uses AlphaFold's AMBER minimization on crystal PDBs

set -eo pipefail

# ACCRE migrated module system — load new stack first
setup_accre_software_stack 2>/dev/null || true
source /sb/apps/alphafold232/miniconda3/bin/activate af232

set -u

DATA_ROOT="/data/p_csb_meiler/agarwm5/af_work"
CLEANED="/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/cleaned"
SCRIPT="/data/p_csb_meiler/agarwm5/protein_pipeline/scripts/amber_relax_crystal_v3.py"

mapfile -t TARGETS < <(ls -1 "$DATA_ROOT" | grep -v "AMBER_FIX" | sort)
IDX=$((SLURM_ARRAY_TASK_ID-1))
TARGET="${TARGETS[$IDX]}"
D="$DATA_ROOT/$TARGET"
CRYSTAL="$CLEANED/${TARGET}.pdb"

echo "AMBER crystal fill: $TARGET on $(hostname)"
nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || true

python "$SCRIPT" "$D" "$CRYSTAL"
