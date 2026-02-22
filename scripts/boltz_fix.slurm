#!/bin/bash --norc
#SBATCH --constraint=csbtmp
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=batch_gpu
#SBATCH --gres=gpu:nvidia_rtx_a6000:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --job-name=boltz_fix
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/boltz_fix_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/boltz_fix_%A_%a.err
#SBATCH --array=1-12

set -euo pipefail

# 12 targets to resubmit:
# 6 header fixes: 1A2K 3BIW 3EO1 3L89 3R9A 3RVW
# 5 empty FASTA regenerated: 4GXU 5HYS BAAD BOYV BP57
# 1 out of range: CP57
TARGETS=(1A2K 3BIW 3EO1 3L89 3R9A 3RVW 4GXU 5HYS BAAD BOYV BP57 CP57)

DATA_ROOT="/data/p_csb_meiler/agarwm5/af_work"
BOLTZ_MINICONDA="/sb/apps/boltz1-v0.4.1/miniconda3"

IDX=$((SLURM_ARRAY_TASK_ID-1))
TARGET="${TARGETS[$IDX]}"
D="$DATA_ROOT/$TARGET"
echo "Processing: $TARGET"

cd "$D"

# Clean old failed output
rm -rf boltz_out/boltz_results_boltz_input

FASTA="boltz_input.fasta"
[[ -f "$FASTA" ]] || { echo "No boltz_input.fasta in $D"; exit 2; }

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
