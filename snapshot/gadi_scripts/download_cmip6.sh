#!/bin/bash
#PBS -N cmip6_download
#PBS -q copyq
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l walltime=04:00:00
#PBS -l wd
#PBS -l storage=scratch/xe2+gdata/xe2
#PBS -P xe2

module load python3

python3 /g/data/xe2/eh4663/scripts/download_cmip6.py