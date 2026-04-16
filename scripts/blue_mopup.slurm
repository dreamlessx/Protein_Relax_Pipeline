#!/bin/bash --norc
#SBATCH --job-name=blue_mopup
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=14-00:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/blue_mopup_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/blue_mopup_%A_%a.err
# SBATCH --array=1-NLINES%100
# ^^^ Uncomment and set NLINES after generating worklist

# Blue mop-up: one array task per (target, model, protocol) combo
# Each task runs only the missing replicates for ONE combo (~5-50 min each)
# Safe to run alongside existing jobs — skip-logic checks before each rep
#
# Setup:
#   1. bash blue_scan_gaps.sh > /data/p_csb_meiler/agarwm5/protein_pipeline/blue_worklist.tsv
#   2. NLINES=$(wc -l < blue_worklist.tsv)
#   3. Edit --array above: --array=1-${NLINES}%100
#   4. sbatch blue_mopup.slurm

set -euo pipefail

WORKLIST="/data/p_csb_meiler/agarwm5/protein_pipeline/blue_worklist.tsv"
RELAX="/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/source/bin/relax.linuxgccrelease"
DB="/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/database"
DATA_ROOT="/data/p_csb_meiler/agarwm5/af_work"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$WORKLIST")
[[ -z "$LINE" ]] && { echo "No work for task ${SLURM_ARRAY_TASK_ID}"; exit 0; }

IFS=$'\t' read -r TARGET PDB SRC_TYPE MODEL_LABEL PROTOCOL MISSING_REPS <<< "$LINE"

echo "=== Blue mop-up task ${SLURM_ARRAY_TASK_ID} ==="
echo "Target: $TARGET, Model: ${SRC_TYPE}_${MODEL_LABEL}, Protocol: $PROTOCOL"
echo "Missing reps: $MISSING_REPS"
echo "Start: $(date)"

OUT_DIR="$DATA_ROOT/$TARGET/rosetta_out/${SRC_TYPE}_${MODEL_LABEL}/${PROTOCOL}"
mkdir -p "$OUT_DIR/log"

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

# Verify input PDB exists
if [[ ! -f "$PDB" ]]; then
  echo "ERROR: Input PDB not found: $PDB"
  exit 1
fi

IFS=',' read -ra REPS <<< "$MISSING_REPS"
DONE=0
for r in "${REPS[@]}"; do
  # Re-check at runtime (another job may have completed this rep)
  if compgen -G "${OUT_DIR}/*_r${r}.pdb.gz" > /dev/null 2>&1; then
    echo "Rep $r already exists — skipping"
    DONE=$((DONE+1))
    continue
  fi

  echo "Running rep $r..."
  "$RELAX" \
    "${COMMON_FLAGS[@]}" \
    -s "$PDB" \
    ${CFG[$PROTOCOL]} \
    -out:path:all "$OUT_DIR" \
    -scorefile relax.fasc \
    -out::suffix "_r${r}" \
    >> "$OUT_DIR/log/r${r}.log" 2>&1 || {
      echo "WARNING: Rosetta failed on rep $r"
      continue
    }
  DONE=$((DONE+1))
done

echo "Completed $DONE/${#REPS[@]} reps for ${SRC_TYPE}_${MODEL_LABEL} $PROTOCOL"
echo "Done: $(date)"
