!/bin/bash
#PBS -N candidate_extract
#PBS -P xe2
#PBS -q normal
#PBS -l walltime=01:00:00
#PBS -l ncpus=8
#PBS -l mem=32gb
#PBS -l storage=gdata/xe2+scratch/xe2
#PBS -j oe

export PATH=/home/272/eh4663/bin:$PATH
module load plink2

cd /g/data/xe2/projects/paneuc2/notes/2026-07-22_pauci-gea/

plink2 \
  --bfile ash_scdp3_sampmiss50_q30_snpmiss40_biallelic_mac10 \
  --allow-extra-chr \
  --extract bed1 /g/data/xe2/eh4663/scripts/candidates_merged.bed \
  --indep-pairwise 10kb 1 0.7 \
  --out /scratch/xe2/eh4663/LFMM/candidates_prune \
  --threads 8

plink2 \
  --bfile ash_scdp3_sampmiss50_q30_snpmiss40_biallelic_mac10 \
  --allow-extra-chr \
  --extract /scratch/xe2/eh4663/LFMM/candidates_prune_very_strict.prune.in \
  --export A \
  --out /scratch/xe2/eh4663/LFMM/candidates_geno \
  --threads 8