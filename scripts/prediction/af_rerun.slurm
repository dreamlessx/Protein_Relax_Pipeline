#!/bin/bash --norc
#SBATCH --constraint=csbtmp
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=batch_gpu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --gres=gpu:nvidia_rtx_a6000:1
#SBATCH --mem=80G
#SBATCH --time=72:00:00
#SBATCH --job-name=af_rerun
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/af_rerun_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/af_rerun_%A_%a.err
#SBATCH --array=1-3

DIRLIST=/data/p_csb_meiler/agarwm5/protein_pipeline/active/af_rerun_targets.txt
AF2_MINICONDA=/sb/apps/alphafold232/miniconda3
AF2_DATADIR=/csbtmp/alphafold-data.230
AF2_REPO=/sb/apps/alphafold232/alphafold

CALCDIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$DIRLIST")
TARGET=$(basename "$CALCDIR")
echo "Re-running: $TARGET"
cd "$CALCDIR" || exit 1

rm -rf $CALCDIR/af_out $CALCDIR/af_out_unrelaxed

FASTA=sequence.fasta
[[ -f "$FASTA" ]] || FASTA=boltz_input.fasta
[[ -f "$FASTA" ]] || { echo "ERROR: No FASTA"; exit 2; }

source $AF2_MINICONDA/bin/activate af232
export LD_LIBRARY_PATH=$AF2_MINICONDA/envs/af232/lib:$LD_LIBRARY_PATH

SEQCOUNT=$(grep -c '^>' "$FASTA")

run_alphafold() {
    local DB_PRESET=$1
    if [[ "$DB_PRESET" == "full_dbs" ]]; then
        DB_FLAGS="--bfd_database_path=$AF2_DATADIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt --uniref30_database_path=$AF2_DATADIR/uniref30/UniRef30_2021_03"
    else
        DB_FLAGS="--small_bfd_database_path=$AF2_DATADIR/small_bfd/bfd-first_non_consensus_sequences.fasta"
    fi

    python $AF2_REPO/run_alphafold.py --fasta_paths=$FASTA --max_template_date=9999-12-31 \
        --data_dir=$AF2_DATADIR --output_dir=$CALCDIR/af_out --db_preset=$DB_PRESET --use_gpu_relax \
        --models_to_relax=all \
        --uniref90_database_path=$AF2_DATADIR/uniref90/uniref90.fasta \
        --mgnify_database_path=$AF2_DATADIR/mgnify/mgy_clusters_2022_05.fa \
        $DB_FLAGS \
        --template_mmcif_dir=$AF2_DATADIR/pdb_mmcif/mmcif_files \
        --pdb_seqres_database_path=$AF2_DATADIR/pdb_seqres/pdb_seqres.txt \
        --obsolete_pdbs_path=$AF2_DATADIR/pdb_mmcif/obsolete.dat \
        --uniprot_database_path=$AF2_DATADIR/uniprot/uniprot.fasta \
        --model_preset=multimer --num_multimer_predictions_per_model=1
    return $?
}

if run_alphafold "full_dbs"; then
    DB_USED="full_dbs"
else
    rm -rf $CALCDIR/af_out
    if run_alphafold "reduced_dbs"; then
        DB_USED="reduced_dbs"
    else
        exit 1
    fi
fi

echo "Database used: $DB_USED" > $CALCDIR/db_preset_used.txt

SEQDIR=$(find $CALCDIR/af_out -maxdepth 1 -type d -name "*" | grep -v "^$CALCDIR/af_out$" | head -1)
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

echo "Done: $TARGET"
