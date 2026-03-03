#!/bin/bash --norc
#SBATCH --job-name=boltz_h100
#SBATCH --account=p_meiler_acc
#SBATCH --partition=batch_gpu
#SBATCH --gres=gpu:nvidia_h100_80gb:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/boltz_h100_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/boltz_h100_%A_%a.err
#SBATCH --array=1-5

set -euo pipefail

# 5 targets that OOMed on A6000 48GB but work on H100 80GB (matching Green)
TARGETS=(1B6C 1F6M 1IB1 1S78 2I9B)

DATA_ROOT="/data/p_csb_meiler/agarwm5/af_work"
BOLTZ_MINICONDA="/sb/apps/boltz1-v0.4.1/miniconda3"

IDX=$((SLURM_ARRAY_TASK_ID-1))
TARGET="${TARGETS[$IDX]}"
D="$DATA_ROOT/$TARGET"
echo "Processing (H100 80GB): $TARGET"

cd "$D"

# Clean old failed output
rm -rf boltz_out/boltz_results_boltz_input

FASTA="boltz_input.fasta"
[[ -s "$FASTA" ]] || FASTA="sequence.fasta"
[[ -s "$FASTA" ]] || { echo "No FASTA in $D"; exit 2; }

# Convert to Boltz format if needed
if ! grep -q '|PROTEIN|' "$FASTA" 2>/dev/null; then
  python3 -c "
chains = []
current_header = ''
current_seq = ''
with open('$FASTA') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if current_seq:
                chains.append((current_header, current_seq))
            current_header = line
            current_seq = ''
        else:
            current_seq += line
    if current_seq:
        chains.append((current_header, current_seq))
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
with open('boltz_input.fasta', 'w') as f:
    for i, (h, seq) in enumerate(chains):
        is_dna = all(c in 'ACGT' for c in seq) and len(seq) < 100
        mol = 'DNA' if is_dna else 'PROTEIN'
        f.write(f'>{letters[i]}|{mol}|\n{seq}\n')
"
  FASTA="boltz_input.fasta"
fi

source "$BOLTZ_MINICONDA/bin/activate" boltz1
export LD_LIBRARY_PATH="$BOLTZ_MINICONDA/envs/boltz1/lib:${LD_LIBRARY_PATH:-}"

mkdir -p boltz_out

echo "Starting Boltz prediction (H100 80GB)..."
boltz predict "$FASTA" \
  --use_msa_server \
  --diffusion_samples 5 \
  --recycling_steps 10 \
  --sampling_steps 200 \
  --output_format pdb \
  --out_dir boltz_out \
  --override

echo "Done: $TARGET"
