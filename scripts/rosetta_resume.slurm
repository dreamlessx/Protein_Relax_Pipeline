#!/bin/bash --norc
#SBATCH --job-name=rosetta_resume
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=120:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/rosetta_resume_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/rosetta_resume_%A_%a.err
#SBATCH --array=4,7,10,12,13,16,17,22,24,28,31,36,37,38,39,41,42

set -euo pipefail

DATA_ROOT="/data/p_csb_meiler/agarwm5/af_work"
RELAX="/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/source/bin/relax.linuxgccrelease"
DB="/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/database"

mapfile -t TARGETS < <(ls -1 "$DATA_ROOT" | sort)
N=${#TARGETS[@]}
IDX=$((SLURM_ARRAY_TASK_ID-1))

if (( IDX < 0 || IDX >= N )); then
  echo "No work for task $SLURM_ARRAY_TASK_ID"
  exit 0
fi

TARGET="${TARGETS[$IDX]}"
D="$DATA_ROOT/$TARGET"
echo "Resuming: $TARGET"

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

shopt -s nullglob
ALL_PDBS=()
SRC_TYPES=()
[[ -d "$D/af_out" ]] && for f in "$D/af_out"/ranked_*.pdb; do ALL_PDBS+=("$f"); SRC_TYPES+=("af_relaxed"); done
[[ -d "$D/af_out_unrelaxed" ]] && for f in "$D/af_out_unrelaxed"/ranked_*.pdb; do ALL_PDBS+=("$f"); SRC_TYPES+=("af_unrelaxed"); done
BOLTZ_DIR="$D/boltz_out/boltz_results_boltz_input/predictions/boltz_input"
[[ -d "$BOLTZ_DIR" ]] && for f in "$BOLTZ_DIR"/boltz_input_model_*.pdb; do ALL_PDBS+=("$f"); SRC_TYPES+=("boltz"); done
shopt -u nullglob

if [[ ${#ALL_PDBS[@]} -eq 0 ]]; then
  echo "No PDBs found for $TARGET"
  exit 0
fi

echo "Found ${#ALL_PDBS[@]} models for $TARGET (resuming incomplete)"

for i in "${!ALL_PDBS[@]}"; do
  pdb="${ALL_PDBS[$i]}"
  src_type="${SRC_TYPES[$i]}"
  model=$(basename "$pdb" .pdb)

  for protocol in cartesian_beta cartesian_ref15 dualspace_beta dualspace_ref15 normal_beta normal_ref15; do
    OUT_DIR="$D/rosetta_out/${src_type}_${model}/${protocol}"

    if [[ -d "$OUT_DIR" ]] && [[ $(ls "$OUT_DIR"/*.pdb.gz 2>/dev/null | wc -l) -ge 5 ]]; then
      continue
    fi

    mkdir -p "$OUT_DIR/log"

    for r in 1 2 3 4 5; do
      [[ -f "$OUT_DIR"/*_r${r}.pdb.gz ]] && continue

      "$RELAX" \
        "${COMMON_FLAGS[@]}" \
        -s "$pdb" \
        ${CFG[$protocol]} \
        -out:path:all "$OUT_DIR" \
        -scorefile relax.fasc \
        -out::suffix "_r${r}" \
        >> "$OUT_DIR/log/r${r}.log" 2>&1 || true
    done
  done
done

echo "Done: $TARGET"
