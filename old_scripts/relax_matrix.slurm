#!/bin/bash --norc
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --job-name=relax-matrix
#SBATCH --output=relax-matrix_%A_%a.out
#SBATCH --array=1-1                 # set to 1-N at submit

set -euo pipefail

ROOT="/dors/meilerlab/home/agarwm5/benchmarking/badset"
RELAX="/dors/meilerlab/apps/rosetta/rosetta-3.14/main/source/bin/relax.linuxgccrelease"
DB="/dors/meilerlab/apps/rosetta/rosetta-3.14/main/database"

# ---- enumerate 4-char subdirs like ABCD, 1a2B, etc.
mapfile -t DIRS < <(find "$ROOT" -mindepth 1 -maxdepth 1 -type d \
  -regextype posix-extended -regex '.*/[0-9A-Za-z]{4}$' | sort)

IDX=$((SLURM_ARRAY_TASK_ID-1))
[[ $IDX -ge 0 && $IDX -lt ${#DIRS[@]} ]] || { echo "No work for task $SLURM_ARRAY_TASK_ID"; exit 0; }

D="${DIRS[$IDX]}"
code="$(basename "$D")"

# ---- choose PDB: prefer ${code}.pdb, else first *.pdb
if [[ -f "$D/$code.pdb" ]]; then
  PDB="$D/$code.pdb"
else
  PDB="$(find "$D" -maxdepth 1 -type f -name '*.pdb' | sort | head -n1)"
fi
[[ -n "${PDB:-}" ]] || { echo "No PDB found in $D"; exit 1; }

echo "Task $SLURM_ARRAY_TASK_ID -> $code : $(basename "$PDB")"

# ---- define 6 protocol/weights combos
# key -> flags (beyond the common flags below)
declare -A CFG
CFG[cart_beta]="-relax:cartesian -beta_nov16 -score:weights beta_nov16_cart"
CFG[cart_ref15]="-relax:cartesian -score:weights ref2015_cart"
CFG[dual_beta]="-relax:dualspace -beta_nov16 -score:weights beta_nov16_cart -nonideal -relax:minimize_bond_angles -relax:minimize_bond_lengths"
CFG[dual_ref15]="-relax:dualspace -score:weights ref2015_cart -nonideal -relax:minimize_bond_angles -relax:minimize_bond_lengths"
CFG[norm_beta]="-beta_nov16 -score:weights beta_nov16"
CFG[norm_ref15]="-score:weights ref2015"

# key -> output folder name
outdir() {
  case "$1" in
    cart_beta)   echo "$D/cartesian_beta" ;;
    cart_ref15)  echo "$D/cartesian_ref15" ;;
    dual_beta)   echo "$D/dualspace_beta" ;;
    dual_ref15)  echo "$D/dualspace_ref15" ;;
    norm_beta)   echo "$D/normal_beta" ;;
    norm_ref15)  echo "$D/normal_ref15" ;;
  esac
}

# ---- common flags (used for all runs)
COMMON_FLAGS=(
  -database "$DB"
  -s "$PDB"
  -ignore_zero_occupancy false
  -nstruct 1
  -no_nstruct_label
  -out:pdb_gz
  -flip_HNQ
  -fa_max_dis 9.0
  -optimization:default_max_cycles 200
  -out:levels all:warning
)

# ---- run 5 replicates per combo
for key in cart_beta cart_ref15 dual_beta dual_ref15 norm_beta norm_ref15; do
  OUT="$(outdir "$key")"
  mkdir -p "$OUT/log"

  echo "==> $code :: $key -> $OUT"
  for r in {1..5}; do
    echo "   - replicate r${r}"
    "$RELAX" \
      "${COMMON_FLAGS[@]}" \
      ${CFG[$key]} \
      -out:path:all "$OUT" \
      -scorefile relax.fasc \
      -out::suffix "_r${r}" \
      >> "$OUT/log/${code}_${key}_r${r}.log" 2>&1
  done
done

echo "Done: $code"
