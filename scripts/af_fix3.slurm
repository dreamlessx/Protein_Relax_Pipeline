#!/bin/bash --norc
#SBATCH --constraint=csbtmp
#SBATCH --account=csb_gpu_acc
#SBATCH --partition=batch_gpu
#SBATCH --gres=gpu:nvidia_rtx_a6000:1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --job-name=af_fix3
#SBATCH --output=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/af_fix3_%A.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/protein_pipeline/logs/af_fix3_%A.err

AF2_MINICONDA=/sb/apps/alphafold232/miniconda3
AF2_DATADIR=/csbtmp/alphafold-data.230
AF2_REPO=/sb/apps/alphafold232/alphafold

TARGET="1FC2"
CALCDIR="/data/p_csb_meiler/agarwm5/af_work/$TARGET"
echo "Re-running: $TARGET (fixed FASTA, cleared old output)"
cd "$CALCDIR"
source $AF2_MINICONDA/bin/activate af232
export LD_LIBRARY_PATH=$AF2_MINICONDA/envs/af232/lib:$LD_LIBRARY_PATH

python $AF2_REPO/run_alphafold.py \
    --fasta_paths=sequence.fasta \
    --max_template_date=2024-01-01 \
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
    --num_multimer_predictions_per_model=1

echo "Done: $TARGET"
