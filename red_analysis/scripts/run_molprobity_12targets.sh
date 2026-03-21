#!/bin/bash
#SBATCH --job-name=red_mp12
#SBATCH --partition=batch
#SBATCH --account=p_csb_meiler
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/data/p_csb_meiler/agarwm5/red_analysis/logs/red_mp12_%j.out
#SBATCH --error=/data/p_csb_meiler/agarwm5/red_analysis/logs/red_mp12_%j.err

PYTHON=/home/agarwm5/miniconda3/envs/red_analysis/bin/python
SCRIPT=/data/p_csb_meiler/agarwm5/red_analysis/scripts/compute_rosetta_molprobity.py

TARGETS="1HIA 1I2M 1I4D 1IB1 1IBR 1IJK 1IRA 1J2J 1JIW 1JK9 1JMO 1JPS"

for T in $TARGETS; do
    echo "$(date): Starting $T"
    $PYTHON -u $SCRIPT blue $T
    echo "$(date): Finished $T (exit code: $?)"
    echo "---"
done

echo "=== ALL 12 TARGETS DONE ==="
echo "Output files:"
for T in $TARGETS; do
    F=/data/p_csb_meiler/agarwm5/red_analysis/metrics/rosetta_molprobity_blue_${T}.tsv
    if [ -f "$F" ]; then
        ROWS=$(wc -l < "$F")
        echo "  $T: $ROWS lines"
    else
        echo "  $T: MISSING"
    fi
done
