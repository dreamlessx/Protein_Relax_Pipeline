#!/bin/bash --norc
#SBATCH --job-name=green_mop2
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=14-00:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/logs/green_mop2_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/logs/green_mop2_%A_%a.err
#SBATCH --array=1-7938%200

set -euo pipefail

WORKLIST="/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/green_worklist_v2.tsv"
RELAX="/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/source/bin/relax.linuxgccrelease"
DB="/data/p_csb_meiler/apps/rosetta/rosetta-3.15/main/database"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$WORKLIST")
[[ -z "$LINE" ]] && exit 0

IFS=$'\t' read -r TARGET PDB SRC_TYPE MODEL_LABEL PROTOCOL MISSING_REPS <<< "$LINE"

DATA_DIR=$(echo "$PDB" | sed 's|/AF/.*||; s|/af_out_unrelaxed/.*||; s|/Boltz/.*||; s|/amber_out/.*||')
[[ "$SRC_TYPE" == "crystal" ]] && DATA_DIR=$(dirname "$PDB")

OUT_DIR="$DATA_DIR/rosetta_out/${SRC_TYPE}_${MODEL_LABEL}/${PROTOCOL}"
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

[[ ! -f "$PDB" ]] && { echo "ERROR: $PDB not found"; exit 1; }

IFS=',' read -ra REPS <<< "$MISSING_REPS"
for r in "${REPS[@]}"; do
  compgen -G "${OUT_DIR}/*_r${r}.pdb.gz" > /dev/null 2>&1 && continue
  "$RELAX" \
    "${COMMON_FLAGS[@]}" \
    -s "$PDB" \
    ${CFG[$PROTOCOL]} \
    -out:path:all "$OUT_DIR" \
    -scorefile relax.fasc \
    -out::suffix "_r${r}" \
    >> "$OUT_DIR/log/r${r}.log" 2>&1 || continue
done
