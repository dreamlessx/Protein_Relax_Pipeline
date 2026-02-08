#!/bin/bash
#
# Array job for Boltz-1 (v0.4.1) structure prediction.
# Runs predictions for all targets in DATA_ROOT, skipping completed ones.
#
# Submit with:  sbatch --array=1-N boltz_array.slurm
#   N = number of targets with FASTA files in DATA_ROOT
#
# Environment overrides:
#   BOLTZ_SAMPLES  - Number of structures per target (default: 5)
#   BOLTZ_RECYCLES - Recycling steps (default: 10)
#   BOLTZ_STEPS    - Sampling steps (default: 200)
#   FORCE=1        - Re-run even if output already exists
#
#SBATCH --job-name=boltz-all
#SBATCH --account=p_meiler_acc
#SBATCH --partition=batch_gpu
#SBATCH --gres=gpu:nvidia_l40s:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --output=boltz_all_%A_%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=1-1

set -euo pipefail

DATA_ROOT="${DATA_ROOT:-/dors/meilerlab/home/agarwm5/benchmarking/data}"
BOLTZ_MINICONDA="/sb/apps/boltz1-v0.4.1/miniconda3"

BOLTZ_SAMPLES="${BOLTZ_SAMPLES:-5}"
BOLTZ_RECYCLES="${BOLTZ_RECYCLES:-10}"
BOLTZ_STEPS="${BOLTZ_STEPS:-200}"
BOLTZ_OUT="${BOLTZ_OUT:-./boltz_out_dir}"
BOLTZ_CACHE="${BOLTZ_CACHE:-$HOME/.boltz}"
FORCE="${FORCE:-0}"

# Build deterministic folder list (prefer boltz_input.fasta)
mapfile -t ALL < <(
  find "$DATA_ROOT" -mindepth 2 -maxdepth 2 -type f \
    \( -name 'boltz_input.fasta' -o -name 'sequence.fasta' \) \
  | sed 's#/[^/]*$##' | sort -u
)

N=${#ALL[@]}
IDX=$((SLURM_ARRAY_TASK_ID-1))
if (( IDX < 0 || IDX >= N )); then
  echo "No work for task $SLURM_ARRAY_TASK_ID (have $N targets)."
  exit 0
fi

D="${ALL[$IDX]}"
FASTA_BASENAME="boltz_input.fasta"
[[ -s "$D/$FASTA_BASENAME" ]] || FASTA_BASENAME="sequence.fasta"
[[ -s "$D/$FASTA_BASENAME" ]] || { echo "No FASTA in $D"; exit 2; }

cd "$D"
echo "Task $SLURM_ARRAY_TASK_ID/$N -> $(pwd)/$FASTA_BASENAME"

# Validate Boltz FASTA headers
if ! awk -F'|' 'BEGIN{bad=0} /^>/{if(NF<2||$2!="PROTEIN"){bad=1; print "Bad header:",$0}} END{exit bad}' "$FASTA_BASENAME"; then
  echo "ERROR: Headers must look like >A|PROTEIN|  (found above)"; exit 3
fi

# Activate env
source "$BOLTZ_MINICONDA/bin/activate" boltz1
boltz --version || true

mkdir -p "$BOLTZ_OUT"

# Skip guard
HAVE=$(find "$BOLTZ_OUT" -path '*/predictions/*/boltz_model_*.pdb' -type f 2>/dev/null | wc -l | awk '{print $1}')
if (( FORCE == 0 && HAVE >= BOLTZ_SAMPLES )); then
  echo "[skip] found $HAVE PDBs (>= $BOLTZ_SAMPLES) in $BOLTZ_OUT"
  exit 0
fi

set -x
boltz predict "$FASTA_BASENAME" \
  --use_msa_server \
  --diffusion_samples "$BOLTZ_SAMPLES" \
  --recycling_steps "$BOLTZ_RECYCLES" \
  --sampling_steps "$BOLTZ_STEPS" \
  --output_format pdb \
  --out_dir "$BOLTZ_OUT" \
  --cache "$BOLTZ_CACHE" \
  --override
set +x

echo "PDBs:"
find "$BOLTZ_OUT" -path '*/predictions/*/boltz_model_*.pdb' -printf "  %P\n" | sort || true
