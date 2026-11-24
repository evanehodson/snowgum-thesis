#!/bin/bash
#PBS -N missing
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=04:00:00
#PBS -l storage=gdata/xe2+scratch/xe2
#PBS -P xe2

export PATH=/home/272/eh4663/bin:$PATH

INPUT="/g/data/xe2/projects/paneuc2/final-variants/2026-05-25_filtered_ash-vim/Eregnans_centurion__shortread-bwa__longread-clair3-glnexus__merged__filtered.vcf.gz"
OUTDIR="/scratch/xe2/eh4663/short-read/1_population_structure/1.2_missing_dist/"

FILENAME=$(basename "$INPUT" .vcf.gz)

plink2 --vcf "$INPUT" \
       --vcf-half-call missing \
       --missing \
       --threads 16 \
       --memory 30000 \
       --out "${OUTDIR}${FILENAME}"