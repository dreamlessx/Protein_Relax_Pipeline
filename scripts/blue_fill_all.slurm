#!/bin/bash --norc
#SBATCH --job-name=blue_fill
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=14-00:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/blue_fill_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/blue_fill_%A_%a.err
#SBATCH --array=1-257%100

# Blue pipeline: fill ALL gaps for one target per task
# Handles all 27 input models (creates dirs as needed)
# 14-day limit, skip-logic for existing output
# Covers: af_relaxed, af_unrelaxed, boltz, amber_af, amber_boltz, crystal, amber_crystal

set -euo pipefail
shopt -s nullglob

DATA_ROOT="/data/p_csb_meiler/agarwm5/af_work"
CRYSTAL_DIR="/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/cleaned"
RELAX="/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/source/bin/relax.linuxgccrelease"
DB="/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/database"

mapfile -t TARGETS < <(ls -1 "$DATA_ROOT" | grep -v "AMBER_FIX" | sort)
IDX=$((SLURM_ARRAY_TASK_ID-1))
[[ $IDX -lt 0 || $IDX -ge ${#TARGETS[@]} ]] && exit 0

TARGET="${TARGETS[$IDX]}"
D="$DATA_ROOT/$TARGET"

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

# Build the full 27-model input list with correct labels
ALL_PDBS=()
SRC_TYPES=()
MODEL_LABELS=()

# 1. AF relaxed (5)
for f in "$D"/af_out/ranked_{0,1,2,3,4}.pdb; do
  [[ -f "$f" ]] || continue
  ALL_PDBS+=("$f"); SRC_TYPES+=("af_relaxed"); MODEL_LABELS+=("$(basename "$f" .pdb)")
done

# 2. AF unrelaxed (5)
for f in "$D"/af_out_unrelaxed/ranked_{0,1,2,3,4}.pdb; do
  [[ -f "$f" ]] || continue
  ALL_PDBS+=("$f"); SRC_TYPES+=("af_unrelaxed"); MODEL_LABELS+=("$(basename "$f" .pdb)")
done

# 3. Boltz (5)
BOLTZ_DIR="$D/boltz_out/boltz_results_boltz_input/predictions/boltz_input"
for f in "$BOLTZ_DIR"/boltz_input_model_{0,1,2,3,4}.pdb; do
  [[ -f "$f" ]] || continue
  ALL_PDBS+=("$f"); SRC_TYPES+=("boltz"); MODEL_LABELS+=("$(basename "$f" .pdb)")
done

# 4. AMBER(AF) with MODEL_LABELS (5)
for i in 0 1 2 3 4; do
  f="$D/amber_out/af_unrelaxed_ranked_${i}/relaxed.pdb"
  [[ -f "$f" ]] && { ALL_PDBS+=("$f"); SRC_TYPES+=("amber_af"); MODEL_LABELS+=("ranked_${i}"); }
done

# 5. AMBER(Boltz) with MODEL_LABELS (5)
for i in 0 1 2 3 4; do
  f="$D/amber_out/boltz_model_${i}/relaxed.pdb"
  [[ -f "$f" ]] && { ALL_PDBS+=("$f"); SRC_TYPES+=("amber_boltz"); MODEL_LABELS+=("model_${i}"); }
done

# 6. Crystal (1)
CRYSTAL_PDB="$CRYSTAL_DIR/${TARGET}.pdb"
[[ -f "$CRYSTAL_PDB" ]] && { ALL_PDBS+=("$CRYSTAL_PDB"); SRC_TYPES+=("crystal"); MODEL_LABELS+=("${TARGET}"); }

# 7. AMBER(Crystal) (1)
for f in "$D"/amber_out/crystal_*/relaxed.pdb "$D"/amber_out/crystal/relaxed.pdb; do
  [[ -f "$f" ]] || continue
  ALL_PDBS+=("$f"); SRC_TYPES+=("amber_crystal"); MODEL_LABELS+=("relaxed")
  break  # only one
done

[[ ${#ALL_PDBS[@]} -eq 0 ]] && { echo "No PDBs for $TARGET"; exit 0; }

# Count work
TOTAL=0; SKIP=0
for i in "${!ALL_PDBS[@]}"; do
  for protocol in cartesian_beta cartesian_ref15 dualspace_beta dualspace_ref15 normal_beta normal_ref15; do
    OUT_DIR="$D/rosetta_out/${SRC_TYPES[$i]}_${MODEL_LABELS[$i]}/${protocol}"
    for r in 1 2 3 4 5; do
      TOTAL=$((TOTAL+1))
      compgen -G "${OUT_DIR}/*_r${r}.pdb.gz" > /dev/null 2>&1 && SKIP=$((SKIP+1))
    done
  done
done

REMAINING=$((TOTAL-SKIP))
if [[ $REMAINING -eq 0 ]]; then
  echo "$TARGET: all $TOTAL runs complete - skipping"
  exit 0
fi

echo "========================================"
echo "FILL Task ${SLURM_ARRAY_TASK_ID}: ${TARGET}"
echo "  ${#ALL_PDBS[@]} inputs, $TOTAL total, $SKIP done, $REMAINING remaining"
echo "  Start: $(date)"
echo "========================================"

for i in "${!ALL_PDBS[@]}"; do
  pdb="${ALL_PDBS[$i]}"
  src_type="${SRC_TYPES[$i]}"
  model="${MODEL_LABELS[$i]}"

  for protocol in cartesian_beta cartesian_ref15 dualspace_beta dualspace_ref15 normal_beta normal_ref15; do
    OUT_DIR="$D/rosetta_out/${src_type}_${model}/${protocol}"

    if [[ -d "$OUT_DIR" ]] && [[ $(ls "$OUT_DIR"/*.pdb.gz 2>/dev/null | wc -l) -ge 5 ]]; then
      continue
    fi

    mkdir -p "$OUT_DIR/log"

    for r in 1 2 3 4 5; do
      compgen -G "${OUT_DIR}/*_r${r}.pdb.gz" > /dev/null 2>&1 && continue

      "$RELAX" \
        "${COMMON_FLAGS[@]}" \
        -s "$pdb" \
        ${CFG[$protocol]} \
        -out:path:all "$OUT_DIR" \
        -scorefile relax.fasc \
        -out::suffix "_r${r}" \
        >> "$OUT_DIR/log/r${r}.log" 2>&1 || {
          echo "WARNING: failed ${src_type}_${model} $protocol r$r"
          continue
        }
    done
  done
done

echo "========================================"
echo "FILL Done: $TARGET ($(date))"
echo "========================================"
