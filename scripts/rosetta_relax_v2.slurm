#!/bin/bash --norc
#SBATCH --job-name=rosetta_v2
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=72:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/rosetta_v2_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/rosetta_v2_%A_%a.err
#SBATCH --array=1-257%50

# Rosetta relaxation for BM5.5 benchmark
# 5 input types x 5 models x 6 protocols x 5 replicates per target
#
# Input types:
#   1. af_relaxed       - AF built-in AMBER relaxed (control)
#   2. af_unrelaxed     - AF unrelaxed (no AMBER)
#   3. boltz            - Boltz-1 unrelaxed (no AMBER)
#   4. amber_af         - Standalone AMBER relaxed AF unrelaxed (test)
#   5. amber_boltz      - Standalone AMBER relaxed Boltz (test)
#
# Rosetta protocols:
#   cartesian_beta, cartesian_ref15, dualspace_beta,
#   dualspace_ref15, normal_beta, normal_ref15
#
# Submit with: sbatch rosetta_relax_v2.slurm

set -euo pipefail

DATA_ROOT="/data/p_csb_meiler/agarwm5/af_work"
RELAX="/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/source/bin/relax.linuxgccrelease"
DB="/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/database"

mapfile -t TARGETS < <(ls -1 "$DATA_ROOT" | grep -v "AMBER_FIX" | sort)
N=${#TARGETS[@]}
IDX=$((SLURM_ARRAY_TASK_ID-1))

if (( IDX < 0 || IDX >= N )); then
  echo "No work for task $SLURM_ARRAY_TASK_ID"
  exit 0
fi

TARGET="${TARGETS[$IDX]}"
D="$DATA_ROOT/$TARGET"
echo "Processing: $TARGET"

# --- 6 Rosetta protocols ---
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

# --- Collect all input PDBs with their source type ---
shopt -s nullglob
ALL_PDBS=()
SRC_TYPES=()

# 1. AF built-in AMBER relaxed
if [[ -d "$D/af_out" ]]; then
  for f in "$D/af_out"/ranked_*.pdb; do
    ALL_PDBS+=("$f")
    SRC_TYPES+=("af_relaxed")
  done
fi

# 2. AF unrelaxed
if [[ -d "$D/af_out_unrelaxed" ]]; then
  for f in "$D/af_out_unrelaxed"/ranked_*.pdb; do
    ALL_PDBS+=("$f")
    SRC_TYPES+=("af_unrelaxed")
  done
fi

# 3. Boltz unrelaxed
BOLTZ_DIR="$D/boltz_out/boltz_results_boltz_input/predictions/boltz_input"
if [[ -d "$BOLTZ_DIR" ]]; then
  for f in "$BOLTZ_DIR"/boltz_input_model_*.pdb; do
    ALL_PDBS+=("$f")
    SRC_TYPES+=("boltz")
  done
fi

# 4. Standalone AMBER of AF unrelaxed
if [[ -d "$D/amber_out" ]]; then
  for i in 0 1 2 3 4; do
    f="$D/amber_out/af_unrelaxed_ranked_${i}/relaxed.pdb"
    if [[ -f "$f" ]]; then
      ALL_PDBS+=("$f")
      SRC_TYPES+=("amber_af")
    fi
  done
fi

# 5. Standalone AMBER of Boltz
if [[ -d "$D/amber_out" ]]; then
  for i in 0 1 2 3 4; do
    f="$D/amber_out/boltz_model_${i}/relaxed.pdb"
    if [[ -f "$f" ]]; then
      ALL_PDBS+=("$f")
      SRC_TYPES+=("amber_boltz")
    fi
  done
fi

shopt -u nullglob

if [[ ${#ALL_PDBS[@]} -eq 0 ]]; then
  echo "No PDBs found for $TARGET"
  exit 0
fi

echo "Found ${#ALL_PDBS[@]} input models for $TARGET"
echo "  af_relaxed:   $(printf '%s\n' "${SRC_TYPES[@]}" | grep -c af_relaxed)"
echo "  af_unrelaxed: $(printf '%s\n' "${SRC_TYPES[@]}" | grep -c af_unrelaxed)"
echo "  boltz:        $(printf '%s\n' "${SRC_TYPES[@]}" | grep -c '^boltz$')"
echo "  amber_af:     $(printf '%s\n' "${SRC_TYPES[@]}" | grep -c amber_af)"
echo "  amber_boltz:  $(printf '%s\n' "${SRC_TYPES[@]}" | grep -c amber_boltz)"

# --- Run relaxation ---
for i in "${!ALL_PDBS[@]}"; do
  pdb="${ALL_PDBS[$i]}"
  src_type="${SRC_TYPES[$i]}"
  model=$(basename "$pdb" .pdb)

  for protocol in cartesian_beta cartesian_ref15 dualspace_beta dualspace_ref15 normal_beta normal_ref15; do
    OUT_DIR="$D/rosetta_out/${src_type}_${model}/${protocol}"

    # Skip if already done (5 replicates)
    if [[ -d "$OUT_DIR" ]] && [[ $(ls "$OUT_DIR"/*.pdb.gz 2>/dev/null | wc -l) -ge 5 ]]; then
      continue
    fi

    mkdir -p "$OUT_DIR/log"

    echo "==> $TARGET :: ${src_type}_${model} :: $protocol"
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
