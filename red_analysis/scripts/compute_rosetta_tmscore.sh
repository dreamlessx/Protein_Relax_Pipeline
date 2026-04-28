#!/bin/bash
###############################################################################
# compute_rosetta_tmscore.sh — Red Analysis Pipeline
#
# Computes TM-score for Rosetta-relaxed structures against crystal references.
# Handles .pdb.gz files (decompresses on the fly).
#
# 6 protocols × 5 reps × N input types = many structures per target.
# Expected per target:
#   Blue: 4 complete input types × 6 protocols × 5 reps = 120 (AMBER bug)
#         + 2 partial AMBER types = ~150 total
#   Green: 6 input types × 6 protocols × 5 reps = 180
#   + crystal: 6 protocols × 5 reps = 30
#   Total per target per pipeline: ~150-210
#
# COMMENT: This is the main analysis. The Rosetta comparison is the
# whole point of the paper. We need to know:
#   1. Which Rosetta protocol produces structures closest to crystal?
#   2. Does Rosetta improve or degrade the input structures?
#   3. Are results consistent between Blue and Green?
#   4. How does protocol choice interact with input type?
###############################################################################

set -euo pipefail

TMSCORE="/data/p_csb_meiler/agarwm5/red_analysis/tmp/TMscore"
BLUE_BASE="/data/p_csb_meiler/agarwm5/af_work"
GREEN_BASE="/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/data"
CRYSTAL_DIR="/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/cleaned"
TARGET_LIST="/data/p_csb_meiler/agarwm5/red_analysis/target_list.txt"
OUTDIR="/data/p_csb_meiler/agarwm5/red_analysis/metrics"

PIPELINE="${1:-blue}"

PROTOCOLS=("cartesian_beta" "cartesian_ref15" "dualspace_beta" "dualspace_ref15" "normal_beta" "normal_ref15")

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    TARGET=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TARGET_LIST")
else
    TARGET="${2:-1A2K}"
fi

if [ -z "$TARGET" ]; then
    echo "ERROR: No target for index ${SLURM_ARRAY_TASK_ID:-N/A}"
    exit 1
fi

CRYSTAL="${CRYSTAL_DIR}/${TARGET}.pdb"
if [ ! -f "$CRYSTAL" ]; then
    echo "ERROR: Crystal not found: $CRYSTAL"
    exit 1
fi

OUTFILE="${OUTDIR}/rosetta_tmscore_${PIPELINE}_${TARGET}.tsv"

if [ ! -f "$OUTFILE" ]; then
    echo -e "target\tpipeline\tsrc_type\tmodel\tprotocol\trep\trmsd\ttmscore\tgdtts\tgdtha\taligned_len\tseq_len" > "$OUTFILE"
fi

# Parse TMscore output (same as compute_tmscore.sh)
run_tmscore_gz() {
    local pdbgz="$1"
    local src_type="$2"
    local model="$3"
    local protocol="$4"
    local rep="$5"

    if [ ! -f "$pdbgz" ]; then
        return
    fi

    local tmpdir=$(mktemp -d)
    local tmppdb="${tmpdir}/input.pdb"
    gunzip -c "$pdbgz" > "$tmppdb" 2>/dev/null || { rm -rf "$tmpdir"; return; }

    local output
    output=$("$TMSCORE" "$tmppdb" "$CRYSTAL" 2>/dev/null) || { rm -rf "$tmpdir"; return; }

    local rmsd=$(echo "$output" | grep "^RMSD of  the common" | awk -F'=' '{print $2}' | awk '{print $1}')
    local tmscore=$(echo "$output" | grep "^TM-score" | head -1 | awk -F'=' '{print $2}' | awk '{print $1}')
    local gdtts=$(echo "$output" | grep "^GDT-TS-score" | awk -F'=' '{print $2}' | awk '{print $1}')
    local gdtha=$(echo "$output" | grep "^GDT-HA-score" | awk -F'=' '{print $2}' | awk '{print $1}')
    local aligned=$(echo "$output" | grep "^Superposition" | awk -F'=' '{print $NF}' | awk '{print $1}')
    local seqlen=$(echo "$output" | grep "by which all scores are normalized" | awk -F'=' '{print $2}' | awk '{print $1}')

    [ -z "$rmsd" ] && rmsd="NA"
    [ -z "$tmscore" ] && tmscore="NA"
    [ -z "$gdtts" ] && gdtts="NA"
    [ -z "$gdtha" ] && gdtha="NA"
    [ -z "$aligned" ] && aligned="NA"
    [ -z "$seqlen" ] && seqlen="NA"

    echo -e "${TARGET}\t${PIPELINE}\t${src_type}\t${model}\t${protocol}\t${rep}\t${rmsd}\t${tmscore}\t${gdtts}\t${gdtha}\t${aligned}\t${seqlen}" >> "$OUTFILE"

    rm -rf "$tmpdir"
}

echo "=== Red Rosetta TM-score: ${TARGET} (${PIPELINE}) ==="

if [ "$PIPELINE" == "blue" ]; then
    ROSETTA_DIR="${BLUE_BASE}/${TARGET}/rosetta_out"
elif [ "$PIPELINE" == "green" ]; then
    ROSETTA_DIR="${GREEN_BASE}/${TARGET}/rosetta_out"
fi

if [ ! -d "$ROSETTA_DIR" ]; then
    echo "WARN: No rosetta_out dir for $TARGET ($PIPELINE)"
    exit 0
fi

# Iterate over all source_type directories
for src_dir in "$ROSETTA_DIR"/*/; do
    src_name=$(basename "$src_dir")
    # Parse src_type and model from directory name
    # Format: {src_type}_{model} e.g., af_relaxed_ranked_0, boltz_boltz_input_model_0
    # We keep the full dir name as the "model" identifier

    for protocol in "${PROTOCOLS[@]}"; do
        protocol_dir="${src_dir}${protocol}"
        [ ! -d "$protocol_dir" ] && continue

        for pdbgz in "$protocol_dir"/*_r[1-5].pdb.gz; do
            [ ! -f "$pdbgz" ] && continue
            # Extract rep number from filename (*_r1.pdb.gz -> 1)
            rep=$(basename "$pdbgz" | grep -oP '_r\K[0-9]+')
            run_tmscore_gz "$pdbgz" "$src_name" "$src_name" "$protocol" "$rep"
        done
    done
done

echo "=== Done: $(wc -l < "$OUTFILE") rows in $OUTFILE ==="
