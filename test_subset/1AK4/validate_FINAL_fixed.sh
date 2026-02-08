#!/bin/bash

# Use SBGrid reduce (works!)
REDUCE="/programs/x86_64-linux/system/sbgrid_bin/reduce"
PHENIX="/programs/x86_64-linux/system/sbgrid_bin/phenix.molprobity"
BASE_DIR="${HOME}/test_subset"
OUTPUT_DIR="${HOME}/molprobity_validation_final"

mkdir -p ${OUTPUT_DIR}/{csvs,logs}
MASTER_PROGRESS="${OUTPUT_DIR}/master_progress.txt"
echo "0" > ${MASTER_PROGRESS}

echo "========================================="
echo "FINAL VALIDATION - FIXED!"
echo "========================================="
echo "Start: $(date)"
echo "========================================="

validate_structure() {
    local PROTEIN=$1
    local CATEGORY=$2
    local SUBCATEGORY=$3
    local PROTOCOL=$4
    local FILE_PATH=$5
    local FILE_NAME=$6
    local PROTEIN_CSV=$7
    
    # Skip if already done
    if [ -f "$PROTEIN_CSV" ]; then
        if grep -q "^${CATEGORY},${SUBCATEGORY},${PROTOCOL},${FILE_NAME}," "$PROTEIN_CSV" 2>/dev/null; then
            return 0
        fi
    fi
    
    [ ! -f "$FILE_PATH" ] && return 1
    
    local WORK_DIR=$(mktemp -d ${OUTPUT_DIR}/temp_XXXXX)
    
    # Decompress
    if [[ $FILE_PATH == *.gz ]]; then
        gunzip -c ${FILE_PATH} > ${WORK_DIR}/input.pdb 2>/dev/null || { rm -rf ${WORK_DIR}; return 1; }
    else
        cp ${FILE_PATH} ${WORK_DIR}/input.pdb || { rm -rf ${WORK_DIR}; return 1; }
    fi
    
    # Add hydrogens (redirect stderr to suppress warnings)
    ${REDUCE} -FLIP -Quiet ${WORK_DIR}/input.pdb > ${WORK_DIR}/input_H.pdb 2>/dev/null || { rm -rf ${WORK_DIR}; return 1; }
    
    # Run Phenix
    cd ${WORK_DIR}
    timeout 300 ${PHENIX} input_H.pdb keep_hydrogens=True > phenix.log 2>&1
    cd - > /dev/null
    
    # Parse and write SINGLE LINE to CSV
    if [ -f "${WORK_DIR}/molprobity.out" ]; then
        local CLASH=$(grep "Clashscore" ${WORK_DIR}/molprobity.out | awk '{print $NF}')
        local RAMA_FAV=$(grep -A 3 "Ramachandran Plot:" ${WORK_DIR}/molprobity.out | grep "Favored" | awk '{print $(NF-1)}')
        local RAMA_ALLOW=$(grep -A 3 "Ramachandran Plot:" ${WORK_DIR}/molprobity.out | grep "Allowed" | awk '{print $(NF-1)}')
        local RAMA_OUT=$(grep -A 3 "Ramachandran Plot:" ${WORK_DIR}/molprobity.out | grep "Outliers" | awk '{print $(NF-1)}')
        local ROT_OUT=$(grep -A 3 "Rotamer:" ${WORK_DIR}/molprobity.out | grep "Outliers" | awk '{print $(NF-1)}')
        local CB_DEV=$(grep "C-beta deviations" ${WORK_DIR}/molprobity.out | awk '{print $NF}')
        local BOND_RMSZ=$(grep "RMS(bonds)" ${WORK_DIR}/molprobity.out | awk '{print $NF}')
        local ANGLE_RMSZ=$(grep "RMS(angles)" ${WORK_DIR}/molprobity.out | awk '{print $NF}')
        local MP_SCORE=$(grep "MolProbity score" ${WORK_DIR}/molprobity.out | awk '{print $NF}')
        
        # SINGLE LINE OUTPUT
        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
            "${CATEGORY}" "${SUBCATEGORY}" "${PROTOCOL}" "${FILE_NAME}" \
            "${CLASH}" "${RAMA_FAV}" "${RAMA_ALLOW}" "${RAMA_OUT}" \
            "${ROT_OUT}" "${CB_DEV}" "${BOND_RMSZ}" "${ANGLE_RMSZ}" "${MP_SCORE}" >> ${PROTEIN_CSV}
    fi
    
    rm -rf ${WORK_DIR}
    echo $(($(cat ${MASTER_PROGRESS}) + 1)) > ${MASTER_PROGRESS}
}

PROTEINS=($(ls -d ${BASE_DIR}/*/ | xargs -n1 basename | sort))
PROTOCOLS=(cartesian_beta cartesian_ref15 dualspace_beta dualspace_ref15 normal_beta normal_ref15)

PROTEIN_NUM=0
for PROTEIN_DIR in ${BASE_DIR}/*/; do
    PROTEIN=$(basename ${PROTEIN_DIR})
    PROTEIN_NUM=$((PROTEIN_NUM + 1))
    
    PROTEIN_CSV="${OUTPUT_DIR}/csvs/${PROTEIN}_validation.csv"
    if [ ! -f "${PROTEIN_CSV}" ]; then
        echo "Category,Subcategory,Protocol,File,Clashscore,Rama_Favored,Rama_Allowed,Rama_Outliers,Rotamer_Outliers,CBeta_Dev,Bond_RMSZ,Angle_RMSZ,MolProbity_Score" > ${PROTEIN_CSV}
    fi
    
    echo "[${PROTEIN_NUM}/20] ${PROTEIN} (Total: $(cat ${MASTER_PROGRESS}))"
    
    # Original PDB
    [ -f "${PROTEIN_DIR}/${PROTEIN}.pdb" ] && \
        validate_structure "$PROTEIN" "Original" "raw" "none" "${PROTEIN_DIR}/${PROTEIN}.pdb" "${PROTEIN}" "${PROTEIN_CSV}"
    
    # AlphaFold raw
    for i in 0 1 2 3 4; do
        [ -f "${PROTEIN_DIR}/AF/ranked_${i}.pdb" ] && \
            validate_structure "$PROTEIN" "AlphaFold" "raw" "none" "${PROTEIN_DIR}/AF/ranked_${i}.pdb" "ranked_${i}" "${PROTEIN_CSV}"
    done
    
    # Boltz raw
    for i in 0 1 2 3 4; do
        [ -f "${PROTEIN_DIR}/Boltz/boltz_input_model_${i}.pdb" ] && \
            validate_structure "$PROTEIN" "Boltz" "raw" "none" "${PROTEIN_DIR}/Boltz/boltz_input_model_${i}.pdb" "boltz_model_${i}" "${PROTEIN_CSV}"
    done
    
    # Original relaxed
    for PROTOCOL in "${PROTOCOLS[@]}"; do
        for r in 1 2 3 4 5; do
            [ -f "${PROTEIN_DIR}/${PROTOCOL}/${PROTEIN}_r${r}.pdb.gz" ] && \
                validate_structure "$PROTEIN" "Original" "relaxed" "${PROTOCOL}" "${PROTEIN_DIR}/${PROTOCOL}/${PROTEIN}_r${r}.pdb.gz" "${PROTEIN}_r${r}" "${PROTEIN_CSV}"
        done
    done
    
    # AF relaxed
    for MODEL in 0 1 2 3 4; do
        for PROTOCOL in "${PROTOCOLS[@]}"; do
            for r in 1 2 3 4 5; do
                [ -f "${PROTEIN_DIR}/relax/AF/ranked_${MODEL}/${PROTOCOL}/ranked_${MODEL}_r${r}.pdb.gz" ] && \
                    validate_structure "$PROTEIN" "AlphaFold" "relaxed_ranked${MODEL}" "${PROTOCOL}" "${PROTEIN_DIR}/relax/AF/ranked_${MODEL}/${PROTOCOL}/ranked_${MODEL}_r${r}.pdb.gz" "ranked_${MODEL}_r${r}" "${PROTEIN_CSV}"
            done
        done
    done
    
    # Boltz relaxed
    for MODEL in 0 1 2 3 4; do
        for PROTOCOL in "${PROTOCOLS[@]}"; do
            for r in 1 2 3 4 5; do
                [ -f "${PROTEIN_DIR}/relax/Boltz/boltz_input_model_${MODEL}/${PROTOCOL}/boltz_input_model_${MODEL}_r${r}.pdb.gz" ] && \
                    validate_structure "$PROTEIN" "Boltz" "relaxed_model${MODEL}" "${PROTOCOL}" "${PROTEIN_DIR}/relax/Boltz/boltz_input_model_${MODEL}/${PROTOCOL}/boltz_input_model_${MODEL}_r${r}.pdb.gz" "boltz_model_${MODEL}_r${r}" "${PROTEIN_CSV}"
            done
        done
    done
done

echo ""
echo "========================================="
echo "COMPLETE! $(date)"
echo "Total: $(cat ${MASTER_PROGRESS}) structures"
echo "========================================="
