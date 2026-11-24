#!/bin/bash
#PBS -N local_score_extract
#PBS -P xe2
#PBS -q normal
#PBS -l walltime=00:30:00
#PBS -l ncpus=4
#PBS -l mem=16gb
#PBS -l jobfs=5gb
#PBS -l storage=gdata/xe2+scratch/xe2
#PBS -j oe

module load python3/3.11.7
cd /g/data/xe2/eh4663/scripts/

python3 local_score.py