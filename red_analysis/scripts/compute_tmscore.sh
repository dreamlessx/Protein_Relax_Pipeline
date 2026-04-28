#!/bin/bash
###############################################################################
# compute_tmscore.sh - Red Analysis Pipeline
#
# Computes TM-score, RMSD, GDT-TS, GDT-HA for all pre-Rosetta structures
# against crystal references. Designed for SLURM array parallelization.
#
# Input types (pre-Rosetta, all 257/257 complete):
#   1. af_relaxed    - AF2 built-in AMBER (ranked_0..4)
#   2. af_unrelaxed  - AF2 unrelaxed (ranked_0..4)
#   3. boltz         - Boltz-1 (model_0..4)
#   4. amber_af      - Standalone AMBER on AF (5 models)
#   5. amber_boltz   - Standalone AMBER on Boltz (5 models)
#
# Output: TSV with columns:
#   target, source, model_idx, rmsd, tmscore, gdtts, gdtha, aligned_len, seq_len
#
# NOTE: Crystal structures are the reference. We compute TM-score normalized
# by the REFERENCE length (crystal), which is the standard for benchmarking
# predicted structures against experimentally determined ones.
###############################################################################

set -euo pipefail

# --- Config ---
TMSCORE="/data/p_csb_meiler/agarwm5/red_analysis/tmp/TMscore"
BLUE_BASE="/data/p_csb_meiler/agarwm5/af_work"
GREEN_BASE="/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/data"
CRYSTAL_DIR="/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/cleaned"
TARGET_LIST="/data/p_csb_meiler/agarwm5/red_analysis/target_list.txt"
OUTDIR="/data/p_csb_meiler/agarwm5/red_analysis/metrics"

# Which pipeline to analyze: "blue" or "green"
PIPELINE="${1:-blue}"

# SLURM array index -> target mapping
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    TARGET=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TARGET_LIST")
else
    # Single target mode for testing
    TARGET="${2:-1A2K}"
fi

if [ -z "$TARGET" ]; then
    echo "ERROR: No target found for index ${SLURM_ARRAY_TASK_ID:-N/A}"
    exit 1
fi

CRYSTAL="${CRYSTAL_DIR}/${TARGET}.pdb"
if [ ! -f "$CRYSTAL" ]; then
    echo "ERROR: Crystal not found: $CRYSTAL"
    exit 1
fi

OUTFILE="${OUTDIR}/tmscore_${PIPELINE}_${TARGET}.tsv"

# Header (only if file doesn't exist)
if [ ! -f "$OUTFILE" ]; then
    echo -e "target\tpipeline\tsource\tmodel_idx\trmsd\ttmscore\tgdtts\tgdtha\taligned_len\tseq_len" > "$OUTFILE"
fi

# --- Function: run TMscore and parse output ---
run_tmscore() {
    local pdb="$1"
    local source="$2"
    local model_idx="$3"

    if [ ! -f "$pdb" ]; then
        echo "WARN: Missing $pdb" >&2
        return
    fi

    # TMscore model.pdb reference.pdb
    # We use crystal as reference (2nd arg) so TM-score is normalized by crystal length
    local output
    output=$("$TMSCORE" "$pdb" "$CRYSTAL" 2>/dev/null)

    # TMscore output format:
    #   RMSD of  the common residues=    0.796
    #   TM-score    = 0.9858  (d0= 4.12)
    #   GDT-TS-score= 0.9940 %(d<1)=...
    #   GDT-HA-score= 0.9657 %(d<0.5)=...
    #   Structure2: ...    Length=  124 (by which all scores are normalized)
    #   Superposition in the TM-score: Length(d<5.0)= 123
    local rmsd=$(echo "$output" | grep "^RMSD of  the common" | awk -F'=' '{print $2}' | awk '{print $1}')
    local tmscore=$(echo "$output" | grep "^TM-score" | head -1 | awk -F'=' '{print $2}' | awk '{print $1}')
    local gdtts=$(echo "$output" | grep "^GDT-TS-score" | awk -F'=' '{print $2}' | awk '{print $1}')
    local gdtha=$(echo "$output" | grep "^GDT-HA-score" | awk -F'=' '{print $2}' | awk '{print $1}')
    local aligned=$(echo "$output" | grep "^Superposition" | awk -F'=' '{print $NF}' | awk '{print $1}')
    local seqlen=$(echo "$output" | grep "by which all scores are normalized" | awk -F'=' '{print $2}' | awk '{print $1}')

    # Fallback if parsing fails
    [ -z "$rmsd" ] && rmsd="NA"
    [ -z "$tmscore" ] && tmscore="NA"
    [ -z "$gdtts" ] && gdtts="NA"
    [ -z "$gdtha" ] && gdtha="NA"
    [ -z "$aligned" ] && aligned="NA"
    [ -z "$seqlen" ] && seqlen="NA"

    echo -e "${TARGET}\t${PIPELINE}\t${source}\t${model_idx}\t${rmsd}\t${tmscore}\t${gdtts}\t${gdtha}\t${aligned}\t${seqlen}" >> "$OUTFILE"
}

# --- Decompress helper for .pdb.gz ---
decompress_and_score() {
    local pdbgz="$1"
    local source="$2"
    local model_idx="$3"

    if [ ! -f "$pdbgz" ]; then
        echo "WARN: Missing $pdbgz" >&2
        return
    fi

    local tmpdir=$(mktemp -d)
    local tmppdb="${tmpdir}/$(basename "$pdbgz" .gz)"
    gunzip -c "$pdbgz" > "$tmppdb"
    run_tmscore "$tmppdb" "$source" "$model_idx"
    rm -rf "$tmpdir"
}

echo "=== Red Analysis: TM-score for ${TARGET} (${PIPELINE}) ==="

if [ "$PIPELINE" == "blue" ]; then
    BASE="${BLUE_BASE}/${TARGET}"

    # 1. AF relaxed (built-in AMBER) - ranked_0..4
    for i in 0 1 2 3 4; do
        run_tmscore "${BASE}/af_out/ranked_${i}.pdb" "af_relaxed" "$i"
    done

    # 2. AF unrelaxed - ranked_0..4
    # Blue stores unrelaxed in af_out_unrelaxed/
    for i in 0 1 2 3 4; do
        pdb="${BASE}/af_out_unrelaxed/ranked_${i}.pdb"
        if [ ! -f "$pdb" ]; then
            # Fallback: some might be in af_out/ with unrelaxed prefix
            pdb="${BASE}/af_out/unrelaxed_model_${i}_ptm.pdb"
        fi
        run_tmscore "$pdb" "af_unrelaxed" "$i"
    done

    # 3. Boltz - model_0..4
    for i in 0 1 2 3 4; do
        run_tmscore "${BASE}/boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_${i}.pdb" "boltz" "$i"
    done

    # 4. Standalone AMBER AF - af_unrelaxed_ranked_0..4/relaxed.pdb
    for i in 0 1 2 3 4; do
        run_tmscore "${BASE}/amber_out/af_unrelaxed_ranked_${i}/relaxed.pdb" "amber_af" "$i"
    done

    # 5. Standalone AMBER Boltz - boltz_model_0..4/relaxed.pdb
    for i in 0 1 2 3 4; do
        run_tmscore "${BASE}/amber_out/boltz_model_${i}/relaxed.pdb" "amber_boltz" "$i"
    done

elif [ "$PIPELINE" == "green" ]; then
    BASE="${GREEN_BASE}/${TARGET}"

    # 1. AF relaxed (built-in AMBER) - AF/ranked_0..4.pdb
    for i in 0 1 2 3 4; do
        run_tmscore "${BASE}/AF/ranked_${i}.pdb" "af_relaxed" "$i"
    done

    # 2. AF unrelaxed - af_out_unrelaxed/ranked_0..4.pdb (symlinks)
    for i in 0 1 2 3 4; do
        run_tmscore "${BASE}/af_out_unrelaxed/ranked_${i}.pdb" "af_unrelaxed" "$i"
    done

    # 3. Boltz - Boltz/boltz_input_model_0..4.pdb
    for i in 0 1 2 3 4; do
        run_tmscore "${BASE}/Boltz/boltz_input_model_${i}.pdb" "boltz" "$i"
    done

    # 4. Standalone AMBER AF - amber_out/af_unrelaxed_0..4/relaxed.pdb
    for i in 0 1 2 3 4; do
        run_tmscore "${BASE}/amber_out/af_unrelaxed_${i}/relaxed.pdb" "amber_af" "$i"
    done

    # 5. Standalone AMBER Boltz - amber_out/boltz_model_0..4/relaxed.pdb
    for i in 0 1 2 3 4; do
        run_tmscore "${BASE}/amber_out/boltz_model_${i}/relaxed.pdb" "amber_boltz" "$i"
    done
fi

echo "=== Done: $OUTFILE ==="
