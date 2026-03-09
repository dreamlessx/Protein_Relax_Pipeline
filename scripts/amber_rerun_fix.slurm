#!/bin/bash --norc
#SBATCH --constraint=csbtmp
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=batch_gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:nvidia_rtx_a6000:1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --job-name=amber_fix
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/amber_fix_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/amber_fix_%A_%a.err

set -uo pipefail

DIRLIST=/data/p_csb_meiler/agarwm5/protein_pipeline/amber_af_unrel_dirlist.txt
SCRIPT=/data/p_csb_meiler/agarwm5/protein_pipeline/scripts/amber_relax.py
AF2_MINICONDA=/sb/apps/alphafold232/miniconda3

CALCDIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$DIRLIST")

if [[ -z "$CALCDIR" ]]; then
    echo "ERROR: Could not read directory for task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

TARGET=$(basename "$CALCDIR")
echo "Task $SLURM_ARRAY_TASK_ID: $TARGET (AMBER rerun fix)"

source $AF2_MINICONDA/bin/activate af232
export LD_LIBRARY_PATH=$AF2_MINICONDA/envs/af232/lib:${LD_LIBRARY_PATH:-}

python "$SCRIPT" "$CALCDIR"

echo "Done: $TARGET"
