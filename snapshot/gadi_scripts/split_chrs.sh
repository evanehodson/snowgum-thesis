#!/bin/bash
#PBS -N split_chrs
#PBS -P xe2
#PBS -q normal
#PBS -l walltime=00:45:00
#PBS -l ncpus=4
#PBS -l mem=16GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/xe2+scratch/xe2

module load bcftools

SRC="/scratch/xe2/eh4663/short-read/0_pre-processing/refilter/Eregnans_centurion__shortread-bwa__longread-clair3-glnexus__merged__refiltered.vcf.gz"
OUT="/scratch/xe2/eh4663/short-read/0_pre-processing/refilter/split_chrs"

bcftools view -r CM074667.1 $SRC -O z -o $OUT/chr_67.vcf.gz && bcftools index $OUT/chr_67.vcf.gz &
bcftools view -r CM074672.1 $SRC -O z -o $OUT/chr_72.vcf.gz && bcftools index $OUT/chr_72.vcf.gz &

wait