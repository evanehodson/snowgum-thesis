#!/bin/bash
#PBS -N plink_candidates
#PBS -P xe2
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -l ncpus=48
#PBS -l mem=192gb
#PBS -l jobfs=10gb
#PBS -l storage=gdata/xe2+scratch/xe2

export PATH=/home/272/eh4663/bin:$PATH

cd /g/data/xe2/projects/paneuc2/notes/2026-07-22_pauci-gea/

plink --bfile ash_scdp3_sampmiss50_q30_snpmiss40_biallelic_mac10 \
  --allow-extra-chr \
  --extract range /scratch/xe2/eh4663/short-read/2_gea/2.1_candidates/candidate_ranges.txt \
  --make-bed \
  --out /scratch/xe2/eh4663/short-read/2_gea/2.1_candidates/candidates
  
wc -l /scratch/xe2/eh4663/short-read/2_gea/2.1_candidates/candidates.bim