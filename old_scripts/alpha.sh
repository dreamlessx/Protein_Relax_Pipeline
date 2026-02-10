#!/bin/bash --norc
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=batch_gpu
#SBATCH --constraint=csbtmp
#SBATCH --mail-user=mudit.agarwal@vanderbilt.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --gres=gpu:nvidia_rtx_a6000:1
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --job-name=af232-array
#SBATCH --output=af232-array_%A_%a.out
#SBATCH --array=1-1   # real range set at submit time

set -euo pipefail
ROOT="/dors/meilerlab/home/agarwm5/benchmarking"

# Build list on the fly: prefer sequence.fasta over boltz_input.fasta
mapfile -t ALL_SEQ < <(find "$ROOT/data" -mindepth 2 -maxdepth 2 -type f -name "sequence.fasta" | sort)
mapfile -t ALL_BOL < <(find "$ROOT/data" -mindepth 2 -maxdepth 2 -type f -name "boltz_input.fasta" | sort)

declare -A pick
for f in "${ALL_BOL[@]}"; do d=$(dirname "$f"); pick["$d"]="$f"; done
for f in "${ALL_SEQ[@]}"; do d=$(dirname "$f"); pick["$d"]="$f"; done  # overwrite with preferred

FILES=("${pick[@]}")
IFS=$'\n' FILES=($(printf "%s\n" "${FILES[@]}" | sort))

IDX=$((SLURM_ARRAY_TASK_ID-1))
[[ $IDX -ge 0 && $IDX -lt ${#FILES[@]} ]] || { echo "No work for task $SLURM_ARRAY_TASK_ID"; exit 0; }

WORKLINE="${FILES[$IDX]}"
WORKDIR="$(dirname "$WORKLINE")"
FASTA_BASENAME="$(basename "$WORKLINE")"

echo "Task $SLURM_ARRAY_TASK_ID -> $WORKDIR/$FASTA_BASENAME"
cd "$WORKDIR"

# Run the folder-local AF script (monomer/multimer aware; multimer set to 5 outputs)
bash "$ROOT/scripts/alphafold_folder.slurm" "$FASTA_BASENAME"
