#!/bin/bash
#===============================================================================
# COMPREHENSIVE MOLPROBITY VALIDATION SCRIPT
# For publication-quality structural validation metrics.
#
# Validates experimental structures, raw AI predictions, and all relaxed
# structures across 6 protocols x 5 replicates. Creates per-protein CSV
# files for resumability.
#
# Usage:
#   bash run_molprobity.sh [base_dir] [output_dir]
#
# Defaults:
#   base_dir   = $HOME/test_subset
#   output_dir = $HOME/molprobity_validation_final
#
# Requirements:
#   - reduce (SBGrid or CCP4)
#   - phenix.molprobity (Phenix suite)
#===============================================================================

# Paths (adjust for your system)
REDUCE="${REDUCE:-/programs/x86_64-linux/system/sbgrid_bin/reduce}"
PHENIX="${PHENIX:-/programs/x86_64-linux/system/sbgrid_bin/phenix.molprobity}"
BASE_DIR="${1:-${HOME}/test_subset}"
OUTPUT_DIR="${2:-${HOME}/molprobity_validation_final}"

mkdir -p "${OUTPUT_DIR}/csvs"
mkdir -p "${OUTPUT_DIR}/logs"

# Progress tracking
MASTER_PROGRESS="${OUTPUT_DIR}/master_progress.txt"
[ ! -f "$MASTER_PROGRESS" ] && echo "0" > "$MASTER_PROGRESS"

# CSV Header
CSV_HEADER="Protein,Category,Subcategory,Protocol,Replicate,Clashscore,Clashscore_Percentile,Ramachandran_Favored,Ramachandran_Allowed,Ramachandran_Outliers,Ramachandran_Outliers_List,Rotamer_Outliers,Rotamer_Outliers_Pct,CBeta_Deviations,CBeta_Outliers_List,Bond_RMSZ,Bond_Outliers,Angle_RMSZ,Angle_Outliers,MolProbity_Score,MolProbity_Percentile"

PROTOCOLS=(cartesian_beta cartesian_ref15 dualspace_beta dualspace_ref15 normal_beta normal_ref15)

echo "========================================="
echo "COMPREHENSIVE MOLPROBITY VALIDATION"
echo "========================================="
echo "Start: $(date)"
echo "Base directory: ${BASE_DIR}"
echo "Output: ${OUTPUT_DIR}"
echo "========================================="

#===============================================================================
# VALIDATION FUNCTION
#===============================================================================
validate_structure() {
    local PROTEIN=$1
    local CATEGORY=$2
    local SUBCATEGORY=$3
    local PROTOCOL=$4
    local REPLICATE=$5
    local FILE_PATH=$6
    local PROTEIN_CSV=$7

    local KEY="${PROTEIN},${CATEGORY},${SUBCATEGORY},${PROTOCOL},${REPLICATE}"

    # Skip if already validated
    if grep -q "^${KEY}," "${PROTEIN_CSV}" 2>/dev/null; then
        return 0
    fi

    if [ ! -f "$FILE_PATH" ]; then
        return 1
    fi

    local WORK_DIR=$(mktemp -d "${OUTPUT_DIR}/temp_XXXXXX")

    # Decompress if needed
    if [[ "$FILE_PATH" == *.gz ]]; then
        gunzip -c "$FILE_PATH" > "${WORK_DIR}/input.pdb" 2>/dev/null
        if [ $? -ne 0 ]; then
            rm -rf "$WORK_DIR"
            return 1
        fi
    else
        cp "$FILE_PATH" "${WORK_DIR}/input.pdb"
    fi

    # Add hydrogens with reduce
    ${REDUCE} -FLIP -Quiet "${WORK_DIR}/input.pdb" > "${WORK_DIR}/input_H.pdb" 2>/dev/null
    if [ ! -s "${WORK_DIR}/input_H.pdb" ]; then
        ${REDUCE} -Quiet "${WORK_DIR}/input.pdb" > "${WORK_DIR}/input_H.pdb" 2>/dev/null
    fi
    if [ ! -s "${WORK_DIR}/input_H.pdb" ]; then
        cp "${WORK_DIR}/input.pdb" "${WORK_DIR}/input_H.pdb"
    fi

    # Run Phenix MolProbity
    cd "${WORK_DIR}"
    timeout 600 ${PHENIX} input_H.pdb keep_hydrogens=True > phenix.log 2>&1
    local PHENIX_EXIT=$?
    cd - > /dev/null

    local OUTFILE="${WORK_DIR}/molprobity.out"

    if [ -f "$OUTFILE" ]; then
        local CLASHSCORE=$(grep "Clashscore" "$OUTFILE" | head -1 | awk '{print $NF}')
        local CLASH_PCTL=$(grep "Clashscore.*percentile" "$OUTFILE" | grep -oP '\d+(?=th|\()' | head -1)
        [ -z "$CLASH_PCTL" ] && CLASH_PCTL="NA"

        local RAMA_FAV=$(grep -A5 "Ramachandran" "$OUTFILE" | grep -i "favored" | awk '{for(i=1;i<=NF;i++) if($i ~ /[0-9]+\.[0-9]+%?/) print $i}' | tr -d '%' | head -1)
        local RAMA_ALLOW=$(grep -A5 "Ramachandran" "$OUTFILE" | grep -i "allowed" | awk '{for(i=1;i<=NF;i++) if($i ~ /[0-9]+\.[0-9]+%?/) print $i}' | tr -d '%' | head -1)
        local RAMA_OUT=$(grep -A5 "Ramachandran" "$OUTFILE" | grep -i "outlier" | awk '{for(i=1;i<=NF;i++) if($i ~ /[0-9]+\.[0-9]+%?/) print $i}' | tr -d '%' | head -1)
        local RAMA_OUT_LIST=$(grep -A50 "Ramachandran outliers" "$OUTFILE" | grep -E "^\s+[A-Z]{3}\s+[A-Z]" | head -5 | tr '\n' ';' | sed 's/;$//')
        [ -z "$RAMA_OUT_LIST" ] && RAMA_OUT_LIST="None"

        local ROT_OUT_NUM=$(grep -A5 "Rotamer" "$OUTFILE" | grep -i "outlier" | awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+$/) print $i}' | head -1)
        local ROT_OUT_PCT=$(grep -A5 "Rotamer" "$OUTFILE" | grep -i "outlier" | awk '{for(i=1;i<=NF;i++) if($i ~ /[0-9]+\.[0-9]+%?/) print $i}' | tr -d '%' | head -1)

        local CBETA_DEV=$(grep -i "C-beta" "$OUTFILE" | awk '{print $NF}')
        local CBETA_LIST=$(grep -A20 "C-beta deviations" "$OUTFILE" | grep -E "^\s+[A-Z]{3}\s+[A-Z]" | head -5 | tr '\n' ';' | sed 's/;$//')
        [ -z "$CBETA_LIST" ] && CBETA_LIST="None"

        local BOND_RMSZ=$(grep -i "RMS.bonds" "$OUTFILE" | awk '{print $NF}')
        local BOND_OUT=$(grep -i "bonds.*outlier" "$OUTFILE" | awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+$/) print $i}' | head -1)
        [ -z "$BOND_OUT" ] && BOND_OUT="0"

        local ANGLE_RMSZ=$(grep -i "RMS.angles" "$OUTFILE" | awk '{print $NF}')
        local ANGLE_OUT=$(grep -i "angles.*outlier" "$OUTFILE" | awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+$/) print $i}' | head -1)
        [ -z "$ANGLE_OUT" ] && ANGLE_OUT="0"

        local MP_SCORE=$(grep -i "MolProbity score" "$OUTFILE" | awk '{print $NF}')
        local MP_PCTL=$(grep -i "MolProbity score.*percentile" "$OUTFILE" | grep -oP '\d+(?=th|\()' | head -1)
        [ -z "$MP_PCTL" ] && MP_PCTL="NA"

        [ -z "$CLASHSCORE" ] && CLASHSCORE="NA"
        [ -z "$RAMA_FAV" ] && RAMA_FAV="NA"
        [ -z "$RAMA_ALLOW" ] && RAMA_ALLOW="NA"
        [ -z "$RAMA_OUT" ] && RAMA_OUT="NA"
        [ -z "$ROT_OUT_NUM" ] && ROT_OUT_NUM="NA"
        [ -z "$ROT_OUT_PCT" ] && ROT_OUT_PCT="NA"
        [ -z "$CBETA_DEV" ] && CBETA_DEV="NA"
        [ -z "$BOND_RMSZ" ] && BOND_RMSZ="NA"
        [ -z "$ANGLE_RMSZ" ] && ANGLE_RMSZ="NA"
        [ -z "$MP_SCORE" ] && MP_SCORE="NA"

        RAMA_OUT_LIST=$(echo "$RAMA_OUT_LIST" | tr ',' ';')
        CBETA_LIST=$(echo "$CBETA_LIST" | tr ',' ';')

        echo "${KEY},${CLASHSCORE},${CLASH_PCTL},${RAMA_FAV},${RAMA_ALLOW},${RAMA_OUT},\"${RAMA_OUT_LIST}\",${ROT_OUT_NUM},${ROT_OUT_PCT},${CBETA_DEV},\"${CBETA_LIST}\",${BOND_RMSZ},${BOND_OUT},${ANGLE_RMSZ},${ANGLE_OUT},${MP_SCORE},${MP_PCTL}" >> "${PROTEIN_CSV}"

        local CURRENT=$(cat "$MASTER_PROGRESS")
        echo $((CURRENT + 1)) > "$MASTER_PROGRESS"
    else
        echo "${KEY},FAILED,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA" >> "${PROTEIN_CSV}"
        echo "FAILED: ${KEY}" >> "${OUTPUT_DIR}/logs/failures.log"
    fi

    rm -rf "$WORK_DIR"
    return 0
}

#===============================================================================
# MAIN LOOP
#===============================================================================
PROTEINS=($(ls -d ${BASE_DIR}/*/ 2>/dev/null | xargs -n1 basename | sort))
TOTAL_PROTEINS=${#PROTEINS[@]}

echo "Found ${TOTAL_PROTEINS} proteins"
echo ""

PROTEIN_NUM=0
for PROTEIN in "${PROTEINS[@]}"; do
    PROTEIN_NUM=$((PROTEIN_NUM + 1))
    PROTEIN_DIR="${BASE_DIR}/${PROTEIN}"
    PROTEIN_CSV="${OUTPUT_DIR}/csvs/${PROTEIN}_validation.csv"

    if [ ! -f "${PROTEIN_CSV}" ]; then
        echo "${CSV_HEADER}" > "${PROTEIN_CSV}"
    fi

    echo "[${PROTEIN_NUM}/${TOTAL_PROTEINS}] Processing ${PROTEIN}..."
    echo "  Progress: $(cat ${MASTER_PROGRESS}) structures validated"

    # 1. Original experimental structure
    if [ -f "${PROTEIN_DIR}/${PROTEIN}.pdb" ]; then
        validate_structure "$PROTEIN" "Experimental" "original" "none" "1" \
            "${PROTEIN_DIR}/${PROTEIN}.pdb" "${PROTEIN_CSV}"
    fi

    # 2. AlphaFold raw predictions (5 models)
    for MODEL in 0 1 2 3 4; do
        if [ -f "${PROTEIN_DIR}/AF/ranked_${MODEL}.pdb" ]; then
            validate_structure "$PROTEIN" "AlphaFold" "raw" "none" "model${MODEL}" \
                "${PROTEIN_DIR}/AF/ranked_${MODEL}.pdb" "${PROTEIN_CSV}"
        fi
    done

    # 3. Boltz-1 raw predictions (5 models)
    for MODEL in 0 1 2 3 4; do
        if [ -f "${PROTEIN_DIR}/Boltz/boltz_input_model_${MODEL}.pdb" ]; then
            validate_structure "$PROTEIN" "Boltz" "raw" "none" "model${MODEL}" \
                "${PROTEIN_DIR}/Boltz/boltz_input_model_${MODEL}.pdb" "${PROTEIN_CSV}"
        fi
    done

    # 4. Experimental structure - relaxed (6 protocols x 5 replicates)
    for PROTOCOL in "${PROTOCOLS[@]}"; do
        for REP in 1 2 3 4 5; do
            if [ -f "${PROTEIN_DIR}/${PROTOCOL}/${PROTEIN}_r${REP}.pdb.gz" ]; then
                validate_structure "$PROTEIN" "Experimental" "relaxed" "$PROTOCOL" "r${REP}" \
                    "${PROTEIN_DIR}/${PROTOCOL}/${PROTEIN}_r${REP}.pdb.gz" "${PROTEIN_CSV}"
            fi
        done
    done

    # 5. AlphaFold - relaxed (5 models x 6 protocols x 5 replicates)
    for MODEL in 0 1 2 3 4; do
        for PROTOCOL in "${PROTOCOLS[@]}"; do
            for REP in 1 2 3 4 5; do
                RELAX_FILE="${PROTEIN_DIR}/relax/AF/ranked_${MODEL}/${PROTOCOL}/ranked_${MODEL}_r${REP}.pdb.gz"
                if [ -f "$RELAX_FILE" ]; then
                    validate_structure "$PROTEIN" "AlphaFold" "relaxed_model${MODEL}" "$PROTOCOL" "r${REP}" \
                        "$RELAX_FILE" "${PROTEIN_CSV}"
                fi
            done
        done
    done

    # 6. Boltz-1 - relaxed (5 models x 6 protocols x 5 replicates)
    for MODEL in 0 1 2 3 4; do
        for PROTOCOL in "${PROTOCOLS[@]}"; do
            for REP in 1 2 3 4 5; do
                RELAX_FILE="${PROTEIN_DIR}/relax/Boltz/boltz_input_model_${MODEL}/${PROTOCOL}/boltz_input_model_${MODEL}_r${REP}.pdb.gz"
                if [ -f "$RELAX_FILE" ]; then
                    validate_structure "$PROTEIN" "Boltz" "relaxed_model${MODEL}" "$PROTOCOL" "r${REP}" \
                        "$RELAX_FILE" "${PROTEIN_CSV}"
                fi
            done
        done
    done

    echo "  Completed ${PROTEIN}: $(wc -l < ${PROTEIN_CSV}) entries"
    echo ""
done

#===============================================================================
# FINAL SUMMARY
#===============================================================================
echo "========================================="
echo "VALIDATION COMPLETE"
echo "========================================="
echo "End time: $(date)"
echo "Total structures validated: $(cat ${MASTER_PROGRESS})"
echo ""
echo "Per-protein CSVs: ${OUTPUT_DIR}/csvs/"
echo "Failure log: ${OUTPUT_DIR}/logs/failures.log"
echo ""
echo "To merge all CSVs into one file:"
echo "  head -1 ${OUTPUT_DIR}/csvs/\$(ls ${OUTPUT_DIR}/csvs/ | head -1) > ${OUTPUT_DIR}/all_validation.csv"
echo "  tail -n +2 -q ${OUTPUT_DIR}/csvs/*_validation.csv >> ${OUTPUT_DIR}/all_validation.csv"
echo "========================================="
