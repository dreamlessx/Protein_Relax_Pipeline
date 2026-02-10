#!/bin/bash --norc
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --job-name=relax-matrix
#SBATCH --output=relax-matrix_%A_%a.out
#SBATCH --array=1-1

set -euo pipefail

# --- paths
ROOT="/dors/meilerlab/home/agarwm5/benchmarking/test"
RELAX="/dors/meilerlab/apps/rosetta/rosetta-3.14/main/source/bin/relax.linuxgccrelease"
DB="/dors/meilerlab/apps/rosetta/rosetta-3.14/main/database"

# --- enumerate 4-char pdb-code dirs (e.g., 1AK4) under ROOT
mapfile -t DIRS < <(find "$ROOT" -mindepth 1 -maxdepth 1 -type d \
  -regextype posix-extended -regex '.*/[0-9A-Za-z]{4}$' | sort)

IDX=$((SLURM_ARRAY_TASK_ID-1))
[[ $IDX -ge 0 && $IDX -lt ${#DIRS[@]} ]] || { echo "No work for task $SLURM_ARRAY_TASK_ID"; exit 0; }

D="${DIRS[$IDX]}"
code="$(basename "$D")"
echo "==> Task $SLURM_ARRAY_TASK_ID on $code"

# ---- define 6 protocol/weights combos
declare -A CFG
CFG[cart_beta]="-relax:cartesian -beta_nov16 -score:weights beta_nov16_cart"
CFG[cart_ref15]="-relax:cartesian -score:weights ref2015_cart"
CFG[dual_beta]="-relax:dualspace -beta_nov16 -score:weights beta_nov16_cart -nonideal -relax:minimize_bond_angles -relax:minimize_bond_lengths"
CFG[dual_ref15]="-relax:dualspace -score:weights ref2015_cart -nonideal -relax:minimize_bond_angles -relax:minimize_bond_lengths"
CFG[norm_beta]="-beta_nov16 -score:weights beta_nov16"
CFG[norm_ref15]="-score:weights ref2015"

# subfolder per mode
outdir_for_key() {
  local base="$1" key="$2"
  case "$key" in
    cart_beta)   echo "$base/cartesian_beta" ;;
    cart_ref15)  echo "$base/cartesian_ref15" ;;
    dual_beta)   echo "$base/dualspace_beta" ;;
    dual_ref15)  echo "$base/dualspace_ref15" ;;
    norm_beta)   echo "$base/normal_beta" ;;
    norm_ref15)  echo "$base/normal_ref15" ;;
  esac
}

# ---- common flags
COMMON_FLAGS=(
  -database "$DB"
  -ignore_zero_occupancy false
  -nstruct 1
  -no_nstruct_label
  -out:pdb_gz
  -flip_HNQ
  -fa_max_dis 9.0
  -optimization:default_max_cycles 200
  -out:levels all:warning
)

# ---- loop over sources (AF, Boltz) and all PDBs within each
for SRC in AF Boltz; do
  [[ -d "$D/$SRC" ]] || { echo "Skipping $SRC (not found in $code)"; continue; }

  mapfile -t PDBS < <(find "$D/$SRC" -maxdepth 1 -type f -name '*.pdb' | sort)
  [[ ${#PDBS[@]} -gt 0 ]] || { echo "No PDBs in $code/$SRC"; continue; }

  for PDB in "${PDBS[@]}"; do
    pdbstem="$(basename "${PDB%.pdb}")"
    # unique base per code / source / pdb name
    OUT_BASE="$D/relax/${SRC}/${pdbstem}"
    mkdir -p "$OUT_BASE/log"

    echo "---- $code :: $SRC :: $(basename "$PDB")"
    for key in cart_beta cart_ref15 dual_beta dual_ref15 norm_beta norm_ref15; do
      OUT="$(outdir_for_key "$OUT_BASE" "$key")"
      mkdir -p "$OUT"

      for r in {1..5}; do
        echo "   [$key] r${r}"
        "$RELAX" \
          -s "$PDB" \
          "${COMMON_FLAGS[@]}" \
          ${CFG[$key]} \
          -out:path:all "$OUT" \
          -scorefile "$OUT_BASE/relax_${pdbstem}.fasc" \
          -out::suffix "_r${r}" \
          >> "$OUT_BASE/log/${code}_${SRC}_${pdbstem}_${key}_r${r}.log" 2>&1
      done
    done
  done
done

echo "Done: $code"
