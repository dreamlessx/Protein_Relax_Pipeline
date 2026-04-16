#!/bin/bash --norc
#SBATCH --job-name=green_axtal_v4
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/green_axtal_v4_%A_%a.out
#SBATCH --error=logs/green_axtal_v4_%A_%a.err
#SBATCH --array=1,7,16,17,21,22,23,25,26,27,28,29,31,41,43,45,46,48,49,50,55,58,59,70,72,73,79,83,87,90,93,101,106,110,111,112,113,116,120,122,123,127,132,136,137,139,141,142,143,146,147,151,155,162,164,165,167,168,176,177,181,182,183,184,185,186,187,188,189,190,192,194,195,196,197,198,199,201,207,209,212,225,239,247,250%100
#SBATCH --exclude=cn1340

# Green Rosetta on AMBER(crystal) outputs

set -euo pipefail

BASE="/panfs/accrepfs.vampire/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking"
DIRLIST="${BASE}/af_dirlist.txt"
RELAX="/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/source/bin/relax.linuxgccrelease"
DB="/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/database"

mkdir -p "${BASE}/logs" 2>/dev/null || true

D=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$DIRLIST")
[[ -z "$D" ]] && { echo "ERROR: No target for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

TARGET=$(basename "$D")

AMBER_CRYSTAL="$D/amber_out/crystal/relaxed.pdb"
if [[ ! -f "$AMBER_CRYSTAL" ]]; then
  echo "No AMBER(crystal) for $TARGET — skipping"
  exit 0
fi

OUT_BASE="$D/rosetta_out/amber_crystal_relaxed"
existing=$(find "$OUT_BASE" -name "*.pdb.gz" 2>/dev/null | wc -l || echo 0)
if [[ "$existing" -ge 30 ]]; then
  echo "Already complete: $TARGET ($existing files)"
  exit 0
fi

echo "Rosetta on AMBER(crystal): $TARGET ($existing/30) ($(date))"

declare -A CFG
CFG[cartesian_beta]="-relax:cartesian -beta_nov16 -score:weights beta_nov16_cart"
CFG[cartesian_ref15]="-relax:cartesian -score:weights ref2015_cart"
CFG[dualspace_beta]="-relax:dualspace -beta_nov16 -score:weights beta_nov16_cart -nonideal -relax:minimize_bond_angles -relax:minimize_bond_lengths"
CFG[dualspace_ref15]="-relax:dualspace -score:weights ref2015_cart -nonideal -relax:minimize_bond_angles -relax:minimize_bond_lengths"
CFG[normal_beta]="-beta_nov16 -score:weights beta_nov16"
CFG[normal_ref15]="-score:weights ref2015"

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

for protocol in cartesian_beta cartesian_ref15 dualspace_beta dualspace_ref15 normal_beta normal_ref15; do
  OUT_DIR="$OUT_BASE/${protocol}"
  done_count=$(ls "$OUT_DIR"/*.pdb.gz 2>/dev/null | wc -l || echo 0)
  if [[ "$done_count" -ge 5 ]]; then continue; fi
  mkdir -p "$OUT_DIR/log"
  for r in 1 2 3 4 5; do
    if compgen -G "${OUT_DIR}/*_r${r}.pdb.gz" > /dev/null 2>&1; then continue; fi
    "$RELAX" "${COMMON_FLAGS[@]}" -s "$AMBER_CRYSTAL" ${CFG[$protocol]} \
      -out:path:all "$OUT_DIR" -scorefile relax.fasc -out::suffix "_r${r}" \
      >> "$OUT_DIR/log/r${r}.log" 2>&1 || {
        echo "WARNING: Rosetta failed on $TARGET amber_crystal, $protocol, rep $r"
        continue
      }
  done
done

echo "Done: $TARGET AMBER(crystal) Rosetta ($(date))"
