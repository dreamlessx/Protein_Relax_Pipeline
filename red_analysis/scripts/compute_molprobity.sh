#!/bin/bash
###############################################################################
# compute_molprobity.sh - Red Analysis Pipeline
#
# Computes MolProbity validation metrics for all pre-Rosetta structures.
# Uses Phenix's molprobity command-line tool.
#
# Metrics extracted:
#   - Clashscore (steric clashes per 1000 atoms)
#   - Ramachandran outliers (%)
#   - Ramachandran favored (%)
#   - Rotamer outliers (%)
#   - MolProbity score (composite)
#   - C-beta deviations
#
# COMMENT: MolProbity is the gold standard for protein structure validation.
# Unlike TM-score (which measures similarity to crystal), MolProbity measures
# intrinsic structural quality regardless of reference. A structure can have
# perfect TM-score but terrible MolProbity if it has bad geometry.
#
# QUESTION: Should we also run on crystal structures themselves? This gives
# a baseline - published crystal structures aren't perfect either. Knowing
# crystal MolProbity helps contextualize prediction quality.
# DECISION: Yes, run on crystals too. It's 257 extra structures, negligible cost.
###############################################################################

set -euo pipefail

# --- Config ---
BLUE_BASE="/data/p_csb_meiler/agarwm5/af_work"
GREEN_BASE="/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/data"
CRYSTAL_DIR="/data/p_csb_meiler/agarwm5/protein_ideal_test/benchmarking/cleaned"
TARGET_LIST="/data/p_csb_meiler/agarwm5/red_analysis/target_list.txt"
OUTDIR="/data/p_csb_meiler/agarwm5/red_analysis/metrics"

PIPELINE="${1:-blue}"

# Source Phenix environment
source /sb/apps/phenix-1.1a/phenix_env 2>/dev/null || {
    echo "ERROR: Cannot source Phenix environment"
    exit 1
}

# SLURM array index -> target mapping
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    TARGET=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TARGET_LIST")
else
    TARGET="${2:-1A2K}"
fi

if [ -z "$TARGET" ]; then
    echo "ERROR: No target for index ${SLURM_ARRAY_TASK_ID:-N/A}"
    exit 1
fi

OUTFILE="${OUTDIR}/molprobity_${PIPELINE}_${TARGET}.tsv"

# Header
if [ ! -f "$OUTFILE" ]; then
    echo -e "target\tpipeline\tsource\tmodel_idx\tclashscore\tramachandran_outliers\tramachandran_favored\trotamer_outliers\tmolprobity_score\tcbeta_deviations" > "$OUTFILE"
fi

# --- Function: run MolProbity and parse output ---
run_molprobity() {
    local pdb="$1"
    local source="$2"
    local model_idx="$3"

    if [ ! -f "$pdb" ]; then
        echo "WARN: Missing $pdb" >&2
        return
    fi

    local tmpdir=$(mktemp -d)

    # Run phenix.molprobity
    # COMMENT: phenix.molprobity writes a bunch of files; use tmpdir to contain them
    cd "$tmpdir"
    local output
    output=$(phenix.molprobity "$pdb" 2>/dev/null) || {
        echo "WARN: MolProbity failed for $pdb" >&2
        rm -rf "$tmpdir"
        return
    }

    # Parse the summary output
    # MolProbity output has lines like:
    #   Clashscore            =   4.52
    #   Ramachandran outliers =   0.00 %
    #   Ramachandran favored  =  97.56 %
    #   Rotamer outliers      =   1.23 %
    #   MolProbity score      =   1.89
    #   C-beta deviations     =      0
    local clashscore=$(echo "$output" | grep "Clashscore" | head -1 | awk -F'=' '{print $2}' | awk '{print $1}')
    local rama_out=$(echo "$output" | grep "Ramachandran outliers" | awk -F'=' '{print $2}' | awk '{print $1}')
    local rama_fav=$(echo "$output" | grep "Ramachandran favored" | awk -F'=' '{print $2}' | awk '{print $1}')
    local rot_out=$(echo "$output" | grep "Rotamer outliers" | awk -F'=' '{print $2}' | awk '{print $1}')
    local mp_score=$(echo "$output" | grep "MolProbity score" | awk -F'=' '{print $2}' | awk '{print $1}')
    local cbeta=$(echo "$output" | grep "C-beta deviations" | awk -F'=' '{print $2}' | awk '{print $1}')

    [ -z "$clashscore" ] && clashscore="NA"
    [ -z "$rama_out" ] && rama_out="NA"
    [ -z "$rama_fav" ] && rama_fav="NA"
    [ -z "$rot_out" ] && rot_out="NA"
    [ -z "$mp_score" ] && mp_score="NA"
    [ -z "$cbeta" ] && cbeta="NA"

    echo -e "${TARGET}\t${PIPELINE}\t${source}\t${model_idx}\t${clashscore}\t${rama_out}\t${rama_fav}\t${rot_out}\t${mp_score}\t${cbeta}" >> "$OUTFILE"

    rm -rf "$tmpdir"
}

echo "=== Red MolProbity: ${TARGET} (${PIPELINE}) ==="

# Crystal structure (baseline)
run_molprobity "${CRYSTAL_DIR}/${TARGET}.pdb" "crystal" "0"

if [ "$PIPELINE" == "blue" ]; then
    BASE="${BLUE_BASE}/${TARGET}"
    for i in 0 1 2 3 4; do
        run_molprobity "${BASE}/af_out/ranked_${i}.pdb" "af_relaxed" "$i"
    done
    for i in 0 1 2 3 4; do
        pdb="${BASE}/af_out_unrelaxed/ranked_${i}.pdb"
        [ ! -f "$pdb" ] && pdb="${BASE}/af_out/unrelaxed_model_${i}_ptm.pdb"
        run_molprobity "$pdb" "af_unrelaxed" "$i"
    done
    for i in 0 1 2 3 4; do
        run_molprobity "${BASE}/boltz_out/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_${i}.pdb" "boltz" "$i"
    done
    for i in 0 1 2 3 4; do
        run_molprobity "${BASE}/amber_out/af_unrelaxed_ranked_${i}/relaxed.pdb" "amber_af" "$i"
    done
    for i in 0 1 2 3 4; do
        run_molprobity "${BASE}/amber_out/boltz_model_${i}/relaxed.pdb" "amber_boltz" "$i"
    done

elif [ "$PIPELINE" == "green" ]; then
    BASE="${GREEN_BASE}/${TARGET}"
    for i in 0 1 2 3 4; do
        run_molprobity "${BASE}/AF/ranked_${i}.pdb" "af_relaxed" "$i"
    done
    for i in 0 1 2 3 4; do
        run_molprobity "${BASE}/af_out_unrelaxed/ranked_${i}.pdb" "af_unrelaxed" "$i"
    done
    for i in 0 1 2 3 4; do
        run_molprobity "${BASE}/Boltz/boltz_input_model_${i}.pdb" "boltz" "$i"
    done
    for i in 0 1 2 3 4; do
        run_molprobity "${BASE}/amber_out/af_unrelaxed_${i}/relaxed.pdb" "amber_af" "$i"
    done
    for i in 0 1 2 3 4; do
        run_molprobity "${BASE}/amber_out/boltz_model_${i}/relaxed.pdb" "amber_boltz" "$i"
    done
fi

echo "=== Done: $OUTFILE ==="
