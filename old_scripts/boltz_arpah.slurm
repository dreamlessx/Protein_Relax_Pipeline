#!/bin/bash --norc
#SBATCH --job-name=boltz-msasrv
#SBATCH --account=p_meiler_acc
#SBATCH --partition=batch_gpu
#SBATCH --gres=gpu:nvidia_l40s:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# --- args: fasta BASENAME in --chdir ---
if [[ $# -ne 1 ]]; then
  echo "Usage: sbatch --chdir=<target_dir> boltz_servermsa.slurm <fasta_basename>" >&2
  exit 2
fi
FASTA="$1"
[[ "$FASTA" == */* ]] && { echo "Provide a basename, not a path: $FASTA" >&2; exit 2; }
[[ -s "$FASTA" ]] || { echo "$FASTA not found in $(pwd)" >&2; exit 2; }

echo "CWD  : $(pwd)"
echo "FASTA: $FASTA"

# --- env ---
BOLTZ1_MINICONDA=/sb/apps/boltz1-v0.4.1/miniconda3
source "$BOLTZ1_MINICONDA/bin/activate" boltz1
boltz --version || true

# --- params (tweakable via env) ---
SAMPLES=${SAMPLES:-5}          # how many structures
RECYCLES=${RECYCLES:-10}
STEPS=${STEPS:-200}
OUTDIR=${OUTDIR:-./boltz_out_dir}
CACHE=${CACHE:-$HOME/.boltz}

mkdir -p "$OUTDIR"

echo "Running: boltz predict (diffusion_samples=$SAMPLES, recycling_steps=$RECYCLES, sampling_steps=$STEPS, output=PDB, use_msa_server)"
set -x
boltz predict "$FASTA" \
  --use_msa_server \
  --diffusion_samples "$SAMPLES" \
  --recycling_steps "$RECYCLES" \
  --sampling_steps "$STEPS" \
  --output_format pdb \
  --out_dir "$OUTDIR" \
  --cache "$CACHE" \
  --override
set +x

echo "PDBs created:"
find "$OUTDIR" -path '*/predictions/*/boltz_model_*.pdb' -printf "  %P\n" | sort || true
COUNT=$(find "$OUTDIR" -path '*/predictions/*/boltz_model_*.pdb' -type f | wc -l | awk '{print $1}')
echo "Total PDBs: $COUNT (expected ${SAMPLES})"
exit 0
