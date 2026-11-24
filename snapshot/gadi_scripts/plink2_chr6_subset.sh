#!/bin/bash
#PBS -N chr6_ld_prep
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l walltime=04:00:00
#PBS -l storage=gdata/xe2+scratch/xe2
#PBS -P xe2

module load bcftools
export PATH=/home/272/eh4663/bin:$PATH

INPUT="/scratch/xe2/eh4663/short-read/0_pre-processing/refilter/Eregnans_centurion__shortread-bwa__longread-clair3-glnexus__merged__refiltered.vcf.gz"
OUTDIR="/scratch/xe2/eh4663/short-read/0_pre-processing/refilter_ch6_subset"
SAMPLES="/g/data/xe2/eh4663/scripts/chr6_all_keep.txt"

mkdir -p "$OUTDIR"

plink2 --vcf "$INPUT" \
  --keep "$SAMPLES" \
  --allow-extra-chr \
  --vcf-half-call missing \
  --chr CM074670.1 --from-bp 38000000 --to-bp 48000000 \
  --thin-count 20000 \
  --threads 16 \
  --export vcf bgz \
  --out "$OUTDIR/chr6_window_downsampled"

bcftools index -c "$OUTDIR/chr6_window_downsampled.vcf.gz"