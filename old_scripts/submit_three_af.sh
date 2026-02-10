#!/usr/bin/env bash
set -euo pipefail

ROOT="/dors/meilerlab/home/agarwm5/benchmarking/data"
SCRIPT="/dors/meilerlab/home/agarwm5/benchmarking/scripts/alphafold_folder.slurm"

# prefer sequence.fasta; fallback to boltz_input.fasta
mapfile -t FILES < <(
  find "$ROOT" -mindepth 2 -maxdepth 2 -type f -name "sequence.fasta" \
  -printf "%h/sequence.fasta\n" 2>/dev/null
)
# if fewer than 3, add boltz_input.fasta
if (( ${#FILES[@]} < 3 )); then
  while IFS= read -r f; do FILES+=("$f"); done < <(
    find "$ROOT" -mindepth 2 -maxdepth 2 -type f -name "boltz_input.fasta" \
    -printf "%h/boltz_input.fasta\n" 2>/dev/null
  )
fi

# unique + first 3
mapfile -t FILES < <(printf "%s\n" "${FILES[@]}" | awk '!seen[$0]++' | sort | head -n 3)

(( ${#FILES[@]} > 0 )) || { echo "No FASTAs found under $ROOT"; exit 2; }

echo "Submitting ${#FILES[@]} AlphaFold jobs..."
for f in "${FILES[@]}"; do
  d=$(dirname "$f")
  id=$(basename "$d")
  base=$(basename "$f")
  echo "  - $id ($base)"
  sbatch --job-name="af2_${id}" --chdir="$d" "$SCRIPT" "$base"
done
