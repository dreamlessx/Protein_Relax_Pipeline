#!/bin/bash --norc
#SBATCH --job-name=boltz_miss
#SBATCH --account=p_meiler_acc
#SBATCH --partition=batch_gpu
#SBATCH --gres=gpu:nvidia_l40s:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=2-00:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/boltz_missing_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/boltz_missing_%A_%a.err

set -euo pipefail

DIRLIST=/data/p_csb_meiler/agarwm5/protein_pipeline/boltz_missing_dirlist.txt
BOLTZ_MINICONDA="/sb/apps/boltz1-v0.4.1/miniconda3"

CALCDIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$DIRLIST")
[[ -z "$CALCDIR" ]] && { echo "ERROR: no dir for task $SLURM_ARRAY_TASK_ID"; exit 1; }

TARGET=$(basename "$CALCDIR")
cd "$CALCDIR"

FASTA="boltz_input.fasta"
[[ -s "$FASTA" ]] || { echo "No boltz_input.fasta in $CALCDIR"; exit 2; }

echo "Task $SLURM_ARRAY_TASK_ID -> $TARGET (missing rerun)"

if [[ -d boltz_out/boltz_results_boltz_input ]]; then
    mv boltz_out/boltz_results_boltz_input boltz_out/boltz_results_old_$(date +%Y%m%d) 2>/dev/null || true
fi

source "$BOLTZ_MINICONDA/bin/activate" boltz1
export LD_LIBRARY_PATH="$BOLTZ_MINICONDA/envs/boltz1/lib:${LD_LIBRARY_PATH:-}"

mkdir -p boltz_out

boltz predict "$FASTA" \
    --use_msa_server \
    --recycling_steps 10 \
    --sampling_steps 200 \
    --diffusion_samples 5 \
    --output_format pdb \
    --out_dir boltz_out \
    --override

echo "Done: $TARGET"
find boltz_out -name "*.pdb" -type f | wc -l
