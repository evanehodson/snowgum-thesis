#!/bin/bash
#PBS -N thin_200k_plink
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=04:00:00
#PBS -l storage=gdata/xe2+scratch/xe2
#PBS -P xe2

export PATH=/home/272/eh4663/bin:$PATH

INPUT="/scratch/xe2/eh4663/short-read/2_gea/2.1_beagle_impute/pauciflora_imputed.vcf.gz"
OUTDIR="/scratch/xe2/eh4663/short-read/0_pre-processing/pauci_beagle_200k_subset/"

mkdir -p "$OUTDIR"

FILENAME=$(basename "$INPUT" .vcf.gz)

plink2 \
  --vcf "$INPUT" \
  --vcf-half-call missing \
  --allow-extra-chr \
  --thin-count 200000 \
  --threads 16 \
  --memory 30000 \
  --export vcf bgz \
  --out "${OUTDIR}${FILENAME}_200k"