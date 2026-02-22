#!/bin/bash --norc
#SBATCH --constraint=csbtmp
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=batch_gpu
#SBATCH --gres=gpu:nvidia_rtx_a6000:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --job-name=boltz_oom
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/boltz_oom_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/boltz_oom_%A_%a.err
#SBATCH --array=1-5

set -euo pipefail

# 5 targets that OOMed for us but worked for Green
# These are in the 1696-2128 residue range
# Try with reduced recycling_steps and diffusion_samples to fit in A6000 48GB
TARGETS=(1B6C 1F6M 1IB1 1S78 2I9B)

DATA_ROOT="/data/p_csb_meiler/agarwm5/af_work"
BOLTZ_MINICONDA="/sb/apps/boltz1-v0.4.1/miniconda3"

IDX=$((SLURM_ARRAY_TASK_ID-1))
TARGET="${TARGETS[$IDX]}"
D="$DATA_ROOT/$TARGET"
echo "Processing (OOM retry): $TARGET"

cd "$D"

# Clean old failed output
rm -rf boltz_out/boltz_results_boltz_input

FASTA="boltz_input.fasta"
[[ -f "$FASTA" ]] || FASTA="sequence.fasta"
[[ -f "$FASTA" ]] || { echo "No FASTA in $D"; exit 2; }

# Convert to Boltz format if needed
if ! grep -q "|PROTEIN|" "$FASTA" 2>/dev/null; then
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

# Reduced settings to avoid OOM on A6000 48GB
boltz predict "$FASTA" \
  --use_msa_server \
  --diffusion_samples 5 \
  --recycling_steps 3 \
  --sampling_steps 200 \
  --output_format pdb \
  --out_dir boltz_out \
  --override

echo "Done: $TARGET"
