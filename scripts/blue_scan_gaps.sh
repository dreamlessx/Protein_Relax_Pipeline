#!/bin/bash
# Blue pipeline: scan for missing Rosetta runs, generate work list
# Output: blue_worklist.tsv (one line per target/model/protocol combo with missing reps)
# Each line: TARGET<tab>PDB_PATH<tab>SRC_TYPE<tab>MODEL_LABEL<tab>PROTOCOL<tab>MISSING_REPS
#
# Usage: bash blue_scan_gaps.sh > blue_worklist.tsv

set -euo pipefail
shopt -s nullglob

DATA_ROOT="/data/p_csb_meiler/agarwm5/af_work"
CRYSTAL_DIR="/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/cleaned"

for D in "$DATA_ROOT"/*/; do
  TARGET=$(basename "$D")
  [[ "$TARGET" == "AMBER_FIX" ]] && continue

  declare -a ALL_PDBS=()
  declare -a SRC_TYPES=()
  declare -a MODEL_LABELS=()

  # AF relaxed
  for f in "$D"af_out/ranked_*.pdb; do
    ALL_PDBS+=("$f")
    SRC_TYPES+=("af_relaxed")
    MODEL_LABELS+=("$(basename "$f" .pdb)")
  done

  # AF unrelaxed
  for f in "$D"af_out_unrelaxed/ranked_*.pdb; do
    ALL_PDBS+=("$f")
    SRC_TYPES+=("af_unrelaxed")
    MODEL_LABELS+=("$(basename "$f" .pdb)")
  done

  # Boltz
  BOLTZ_DIR="${D}boltz_out/boltz_results_boltz_input/predictions/boltz_input"
  for f in "$BOLTZ_DIR"/boltz_input_model_*.pdb; do
    ALL_PDBS+=("$f")
    SRC_TYPES+=("boltz")
    MODEL_LABELS+=("$(basename "$f" .pdb)")
  done

  # AMBER AF (MODEL_LABELS fix)
  for i in 0 1 2 3 4; do
    f="${D}amber_out/af_unrelaxed_ranked_${i}/relaxed.pdb"
    if [[ -f "$f" ]]; then
      ALL_PDBS+=("$f")
      SRC_TYPES+=("amber_af")
      MODEL_LABELS+=("ranked_${i}")
    fi
  done

  # AMBER Boltz (MODEL_LABELS fix)
  for i in 0 1 2 3 4; do
    f="${D}amber_out/boltz_model_${i}/relaxed.pdb"
    if [[ -f "$f" ]]; then
      ALL_PDBS+=("$f")
      SRC_TYPES+=("amber_boltz")
      MODEL_LABELS+=("model_${i}")
    fi
  done

  # Crystal
  CRYSTAL_PDB="$CRYSTAL_DIR/${TARGET}.pdb"
  if [[ -f "$CRYSTAL_PDB" ]]; then
    ALL_PDBS+=("$CRYSTAL_PDB")
    SRC_TYPES+=("crystal")
    MODEL_LABELS+=("${TARGET}")
  fi

  for i in "${!ALL_PDBS[@]}"; do
    pdb="${ALL_PDBS[$i]}"
    src_type="${SRC_TYPES[$i]}"
    model="${MODEL_LABELS[$i]}"

    for protocol in cartesian_beta cartesian_ref15 dualspace_beta dualspace_ref15 normal_beta normal_ref15; do
      OUT_DIR="${D}rosetta_out/${src_type}_${model}/${protocol}"
      missing=""
      for r in 1 2 3 4 5; do
        if ! compgen -G "${OUT_DIR}/*_r${r}.pdb.gz" > /dev/null 2>&1; then
          missing="${missing}${missing:+,}${r}"
        fi
      done
      if [[ -n "$missing" ]]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$TARGET" "$pdb" "$src_type" "$model" "$protocol" "$missing"
      fi
    done
  done

  unset ALL_PDBS SRC_TYPES MODEL_LABELS
done
