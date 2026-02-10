#!/bin/sh
set -euo pipefail

ROOT="/dors/meilerlab/home/agarwm5/benchmarking/data"
SCRIPT="/dors/meilerlab/home/agarwm5/benchmarking/scripts/boltz_arpah.slurm"
INPUT_NAME="boltz_input.fasta"   # use "sequence.fasta" if you skipped step 1




SAMPLES="${SAMPLES:-5}"

# pick any 3 boltz_input.fasta (fall back to sequence.fasta if needed)
mapfile -t FILES < <(find "$ROOT" -mindepth 2 -maxdepth 2 -type f -name "boltz_input.fasta" -printf "%h/boltz_input.fasta\n" | sort | head -n 3)
if (( ${#FILES[@]} < 3 )); then
  while IFS= read -r f; do FILES+=("$f"); done < <(
    find "$ROOT" -mindepth 2 -maxdepth 2 -type f -name "sequence.fasta" -printf "%h/sequence.fasta\n" | sort
  )
  mapfile -t FILES < <(printf "%s\n" "${FILES[@]}" | awk '!seen[$0]++' | head -n 3)
fi

(( ${#FILES[@]} > 0 )) || { echo "No FASTAs found under $ROOT"; exit 2; }
[[ -x "$SCRIPT" ]] || { echo "ERROR: $SCRIPT missing or not executable"; exit 3; }

echo "Submitting ${#FILES[@]} jobs..."
for f in "${FILES[@]}"; do
  d=$(dirname "$f"); id=$(basename "$d"); base=$(basename "$f")
  echo "  - $id ($base)"
  sbatch -A p_meiler_acc -p batch_gpu -t 2-00:00:00 \
    --export=ALL,SAMPLES="$SAMPLES" \
    --job-name="boltz_${id}" --chdir="$d" \
    "$SCRIPT" "$base"
done
