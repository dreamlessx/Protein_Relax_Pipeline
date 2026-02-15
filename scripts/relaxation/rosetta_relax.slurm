#!/bin/bash --norc
#SBATCH --job-name=rosetta_relax
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=72:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/rosetta_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/rosetta_%A_%a.err
#SBATCH --array=1-257%50

set -euo pipefail

DATA_ROOT="/data/p_csb_meiler/agarwm5/af_work"
RELAX="/dors/meilerlab/apps/rosetta/rosetta-3.14/main/source/bin/relax.linuxgccrelease"
DB="/dors/meilerlab/apps/rosetta/rosetta-3.14/main/database"

mapfile -t TARGETS < <(ls -1 "$DATA_ROOT" | sort)
N=${#TARGETS[@]}
IDX=$((SLURM_ARRAY_TASK_ID-1))

if (( IDX < 0 || IDX >= N )); then
  echo "No work for task $SLURM_ARRAY_TASK_ID"
  exit 0
fi

TARGET="${TARGETS[$IDX]}"
D="$DATA_ROOT/$TARGET"
echo "Processing: $TARGET"

# 6 protocols
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

# Collect all AF PDBs (relaxed from af_out, unrelaxed from af_out_unrelaxed)
shopt -s nullglob
ALL_PDBS=()
[[ -d "$D/af_out" ]] && ALL_PDBS+=("$D/af_out"/ranked_*.pdb)
[[ -d "$D/af_out_unrelaxed" ]] && ALL_PDBS+=("$D/af_out_unrelaxed"/ranked_*.pdb)
shopt -u nullglob

if [[ ${#ALL_PDBS[@]} -eq 0 ]]; then
  echo "No AF PDBs found for $TARGET"
  exit 0
fi

echo "Found ${#ALL_PDBS[@]} AF models for $TARGET"

for pdb in "${ALL_PDBS[@]}"; do
  # Determine source type from path
  if [[ "$pdb" == *"/af_out_unrelaxed/"* ]]; then
    src_type="unrelaxed"
  else
    src_type="relaxed"
  fi
  
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
