#!/bin/bash --norc
#SBATCH --constraint=csbtmp
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=batch_gpu
#SBATCH --gres=gpu:nvidia_rtx_a6000:1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --job-name=af_amber_retry
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/af_amber_%A_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/af_amber_%A_%a.err
#SBATCH --array=1-7

TARGETS=(1ATN 1DFJ 1FC2 1WEJ 2BTF 4CPA 5JMO)
AF2_MINICONDA=/sb/apps/alphafold232/miniconda3
AF2_DATADIR=/csbtmp/alphafold-data.230
AF2_REPO=/sb/apps/alphafold232/alphafold

TARGET="${TARGETS[$((SLURM_ARRAY_TASK_ID-1))]}"
CALCDIR="/data/p_csb_meiler/agarwm5/af_work/$TARGET"

echo "Re-running AMBER relaxation for: $TARGET"
cd "$CALCDIR" || exit 1

# Clean previous attempt
rm -rf "$CALCDIR/af_out/sequence/relaxed_*.pdb" "$CALCDIR/af_out/ranked_*.pdb"

FASTA="sequence.fasta"
[[ -f "$FASTA" ]] || FASTA="boltz_input.fasta"
[[ -f "$FASTA" ]] || { echo "ERROR: No FASTA"; exit 2; }

source $AF2_MINICONDA/bin/activate af232
export LD_LIBRARY_PATH=$AF2_MINICONDA/envs/af232/lib:$LD_LIBRARY_PATH

# Run with GPU relax and full error output
python $AF2_REPO/run_alphafold.py \
    --fasta_paths=$FASTA \
    --max_template_date=9999-12-31 \
    --data_dir=$AF2_DATADIR \
    --output_dir=$CALCDIR/af_out \
    --db_preset=reduced_dbs \
    --use_gpu_relax \
    --models_to_relax=all \
    --small_bfd_database_path=$AF2_DATADIR/small_bfd/bfd-first_non_consensus_sequences.fasta \
    --uniref90_database_path=$AF2_DATADIR/uniref90/uniref90.fasta \
    --mgnify_database_path=$AF2_DATADIR/mgnify/mgy_clusters_2022_05.fa \
    --template_mmcif_dir=$AF2_DATADIR/pdb_mmcif/mmcif_files \
    --pdb_seqres_database_path=$AF2_DATADIR/pdb_seqres/pdb_seqres.txt \
    --obsolete_pdbs_path=$AF2_DATADIR/pdb_mmcif/obsolete.dat \
    --uniprot_database_path=$AF2_DATADIR/uniprot/uniprot.fasta \
    --model_preset=multimer \
    --num_multimer_predictions_per_model=1 \
    --use_precomputed_msas=true

# Copy ranked files
SEQDIR="$CALCDIR/af_out/sequence"
if [[ -d "$SEQDIR" ]]; then
    for f in "$SEQDIR"/relaxed_model_*.pdb; do
        [[ -f "$f" ]] || continue
        idx=$(echo "$f" | grep -oP 'model_\K[0-9]+')
        cp "$f" "$CALCDIR/af_out/ranked_$((idx-1)).pdb"
    done
fi

echo "Done: $TARGET"
ls "$CALCDIR/af_out/"
