#!/bin/bash
#SBATCH --job-name=red_initial_energy
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --array=1-257
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/red_analysis/logs/initial_energy_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/red_analysis/logs/initial_energy_%a.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate red_analysis

python /data/p_csb_meiler/agarwm5/red_analysis/scripts/score_initial_energy.py
