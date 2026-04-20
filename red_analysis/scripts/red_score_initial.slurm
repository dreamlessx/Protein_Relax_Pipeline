#!/bin/bash
#SBATCH --job-name=red_score_initial
#SBATCH --account=p_csb_meiler
#SBATCH --partition=batch
#SBATCH --array=1-257
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/red_analysis/logs/score_initial_%a.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/red_analysis/logs/score_initial_%a.err

# Score input PDBs (pre-Rosetta) with Rosetta score_jd2 to get initial energies.
# One target per array task.  Output: metrics/initial_energy/initial_energy_{TARGET}.tsv

source ~/miniconda3/etc/profile.d/conda.sh
conda activate red_analysis

mkdir -p /data/p_csb_meiler/agarwm5/red_analysis/logs
mkdir -p /data/p_csb_meiler/agarwm5/red_analysis/metrics/initial_energy

python /data/p_csb_meiler/agarwm5/red_analysis/scripts/score_initial_energy.py
