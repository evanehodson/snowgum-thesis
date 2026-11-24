#!/bin/bash
#PBS -N merge
#PBS -P xe2
#PBS -l walltime=24:00:00
#PBS -l ncpus=12
#PBS -l mem=64gb
#PBS -l storage=gdata/xe2+scratch/xe2
#PBS -j oe

module load bcftools

mkdir -p /g/data/xe2/projects/paneuc2/final-variants/merged

bcftools merge \
  --write-index \
  -Oz \
  --threads 12 \
  -o /g/data/xe2/projects/paneuc2/final-variants/merged/Eregnans_centurion__shortread-bwa__longread-clair3-glnexus__merged.vcf.gz \
  /g/data/xe2/projects/paneuc2/final-variants/2025-03-16_ash-vim/bwa~Eregnans_centurion~mm_ash.vcf.gz \
  /g/data/xe2/projects/clair3-workflow/final/2026-05-12_clair3-calls/clair3~minimap2~Eregnans_centurion~group_called_glnexus.tags.vcf.gz