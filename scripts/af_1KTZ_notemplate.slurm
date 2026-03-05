#!/bin/bash --norc
#SBATCH --constraint=csbtmp
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=batch_gpu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --gres=gpu:nvidia_rtx_a6000:1
#SBATCH --mem=80G
#SBATCH --time=72:00:00
#SBATCH --job-name=af_1KTZ
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/af_1KTZ_%j.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/af_1KTZ_%j.err

set -uo pipefail

AF2_MINICONDA=/sb/apps/alphafold232/miniconda3
AF2_DATADIR=/csbtmp/alphafold-data.230
AF2_REPO=/sb/apps/alphafold232/alphafold
CALCDIR=/data/p_csb_meiler/agarwm5/af_work/1KTZ
FASTA=sequence.fasta

cd "$CALCDIR" || exit 1
echo "1KTZ: AF re-run with max_template_date=1900-01-01 (skip templates)"

rm -rf af_out 2>/dev/null

source $AF2_MINICONDA/bin/activate af232
export LD_LIBRARY_PATH=$AF2_MINICONDA/envs/af232/lib:${LD_LIBRARY_PATH:-}

SEQCOUNT=$(grep -c '^>' "$FASTA")

run_af() {
    local DB_MODE="$1"
    local COMMON_ARGS=(
        --fasta_paths=$FASTA
        --max_template_date=1900-01-01
        --data_dir=$AF2_DATADIR
        --output_dir=$CALCDIR/af_out
        --use_gpu_relax
        --models_to_relax=all
        --uniref90_database_path=$AF2_DATADIR/uniref90/uniref90.fasta
        --mgnify_database_path=$AF2_DATADIR/mgnify/mgy_clusters_2022_05.fa
        --template_mmcif_dir=$AF2_DATADIR/pdb_mmcif/mmcif_files
        --obsolete_pdbs_path=$AF2_DATADIR/pdb_mmcif/obsolete.dat
    )
    local DB_ARGS=()
    if [[ "$DB_MODE" == "full_dbs" ]]; then
        DB_ARGS=(
            --uniref30_database_path=$AF2_DATADIR/uniref30/UniRef30_2021_03
            --bfd_database_path=$AF2_DATADIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
        )
    else
        DB_ARGS=(
            --db_preset=reduced_dbs
            --small_bfd_database_path=$AF2_DATADIR/small_bfd/bfd-first_non_consensus_sequences.fasta
        )
    fi

    python $AF2_REPO/run_alphafold.py \
        "${COMMON_ARGS[@]}" \
        "${DB_ARGS[@]}" \
        --model_preset=multimer \
        --pdb_seqres_database_path=$AF2_DATADIR/pdb_seqres/pdb_seqres.txt \
        --uniprot_database_path=$AF2_DATADIR/uniprot/uniprot.fasta \
        --num_multimer_predictions_per_model=1
}

if run_af "full_dbs"; then
    echo "full_dbs succeeded"
else
    echo "full_dbs failed, trying reduced_dbs..."
    rm -rf "$CALCDIR/af_out" 2>/dev/null
    run_af "reduced_dbs" || { echo "ERROR: Both failed"; exit 1; }
fi

# Extract ranked PDBs
SEQDIR=$(find $CALCDIR/af_out -maxdepth 1 -type d ! -name af_out | head -1)
if [[ -n "$SEQDIR" ]]; then
    mkdir -p $CALCDIR/af_out_unrelaxed
    for i in 1 2 3 4 5; do
        src=$(ls $SEQDIR/unrelaxed_model_${i}_*.pdb 2>/dev/null | head -1)
        [[ -f "$src" ]] && cp "$src" "$CALCDIR/af_out_unrelaxed/ranked_$((i-1)).pdb"
    done
    for f in $SEQDIR/ranked_*.pdb; do
        [[ -f "$f" ]] && mv "$f" $CALCDIR/af_out/
    done
fi

for subdir in $CALCDIR/af_out/*/; do
    rm -rf "$subdir/msas" 2>/dev/null
    rm -f "$subdir"/features.pkl "$subdir"/result_model_*.pkl "$subdir"/relaxed_model_*.pdb "$subdir"/timings.json 2>/dev/null
done

echo "Done: 1KTZ"
ls $CALCDIR/af_out/ranked_*.pdb 2>/dev/null | wc -l
