#!/bin/bash --norc
#SBATCH --constraint=csbtmp
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=batch_gpu
#SBATCH --gres=gpu:nvidia_rtx_a6000:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --job-name=boltz_bm55
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/boltz_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/boltz_%A_%a.err
#SBATCH --array=1-257%5

set -euo pipefail

DATA_ROOT="/data/p_csb_meiler/agarwm5/af_work"
BOLTZ_MINICONDA="/sb/apps/boltz1-v0.4.1/miniconda3"

mapfile -t TARGETS < <(ls -1 "$DATA_ROOT" | sort)
N=${#TARGETS[@]}
IDX=$((SLURM_ARRAY_TASK_ID-1))

if (( IDX < 0 || IDX >= N )); then
  echo "No work for task $SLURM_ARRAY_TASK_ID"
  exit 0
fi

TARGET="${TARGETS[$IDX]}"
D="$DATA_ROOT/$TARGET"
echo "Processing: $TARGET"

# Check for existing output
if [ -d "$D/boltz_out" ] && [ "$(find "$D/boltz_out" -name "*.pdb" 2>/dev/null | wc -l)" -ge 5 ]; then
  echo "Already complete: $TARGET"
  exit 0
fi

cd "$D"

# Find FASTA
FASTA="boltz_input.fasta"
[[ -f "$FASTA" ]] || FASTA="sequence.fasta"
[[ -f "$FASTA" ]] || { echo "No FASTA in $D"; exit 2; }

# Convert to Boltz format if needed
if ! grep -q "|PROTEIN|" "$FASTA" 2>/dev/null; then
  echo "Converting FASTA to Boltz format..."
  awk '/^>/{
    split($0, a, "|")
    chain = substr(a[1], 2, 1)
    if (chain == "") chain = "A"
    print ">" chain "|PROTEIN|"
    next
  }
  {print}' "$FASTA" > boltz_input.fasta
  FASTA="boltz_input.fasta"
fi

source "$BOLTZ_MINICONDA/bin/activate" boltz1

mkdir -p boltz_out

boltz predict "$FASTA" \
  --use_msa_server \
  --diffusion_samples 5 \
  --recycling_steps 10 \
  --sampling_steps 200 \
  --output_format pdb \
  --out_dir boltz_out \
  --override

echo "Done: $TARGET"
