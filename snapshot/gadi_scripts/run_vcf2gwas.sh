#!/bin/bash
#PBS -N vcf2gwas
#PBS -l ncpus=48
#PBS -l mem=190GB
#PBS -l walltime=24:00:00
#PBS -l storage=gdata/xe2+scratch/xe2
#PBS -P xe2
#PBS -o /scratch/xe2/eh4663/logs/vcf2gwas.out
#PBS -e /scratch/xe2/eh4663/logs/vcf2gwas.err

set -euo pipefail

source /scratch/xe2/eh4663/miniforge3/etc/profile.d/conda.sh
conda activate /g/data/xe2/kevin/conda/envs/paneuc-acanthophis/
source /scratch/xe2/eh4663/cargo/env

cd /scratch/xe2/eh4663/vcf2gwas.snk

snakemake --snakefile vcf2gwas.snk --cores 48 --configfile config.yml