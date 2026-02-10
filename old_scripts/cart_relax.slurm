#!/bin/bash --norc
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00
#SBATCH --job-name=relax-subset
#SBATCH --output=relax-subset_%A_%a.out
#SBATCH --array=1-1

set -euo pipefail

ROOT="/dors/meilerlab/home/agarwm5/benchmarking/subset"
RELAX="/dors/meilerlab/apps/rosetta/rosetta-3.14/main/source/bin/relax.linuxgccrelease"
DB="/dors/meilerlab/apps/rosetta/rosetta-3.14/main/database"

# >>> FIXED FIND <<<
mapfile -t DIRS < <(find "$ROOT" -mindepth 1 -maxdepth 1 -type d \
  -regextype posix-extended -regex '.*/[0-9A-Za-z]{4}$' | sort)

IDX=$((SLURM_ARRAY_TASK_ID-1))
[[ $IDX -ge 0 && $IDX -lt ${#DIRS[@]} ]] || { echo "No work for task $SLURM_ARRAY_TASK_ID"; exit 0; }

D="${DIRS[$IDX]}"
code="$(basename "$D")"

if [[ -f "$D/$code.pdb" ]]; then
  PDB="$D/$code.pdb"
else
  PDB="$(find "$D" -maxdepth 1 -type f -name '*.pdb' | sort | head -n1)"
fi
[[ -n "${PDB:-}" ]] || { echo "No PDB found in $D"; exit 1; }

OUT="$D/cartesian_beta"
mkdir -p "$OUT"

echo "Task $SLURM_ARRAY_TASK_ID -> $code : $(basename "$PDB")"
echo "Output -> $OUT"

"$RELAX" -database "$DB" \
  -s "$PDB" \
  -ignore_zero_occupancy false \
  -nstruct 1 -no_nstruct_label \
  -out:pdb_gz \
  -out:path:all "$OUT" \
  -scorefile relax.fasc \
  -relax:cartesian -beta_nov16 -score:weights beta_nov16_cart \
  -flip_HNQ -fa_max_dis 9.0 \
  -optimization:default_max_cycles 200
