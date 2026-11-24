#!/bin/bash
#PBS -N lfmm_pauci_k9
#PBS -P xe2
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -l ncpus=48
#PBS -l mem=192gb
#PBS -l jobfs=10gb
#PBS -l storage=gdata/xe2+scratch/xe2
#PBS -m abe
#PBS -j oe

set -e

cd /g/data/xe2/projects/paneuc2/notes/2026-07-22_pauci-gea/

bin/schnellfmm \
    -k 9 \
    --bed ash_scdp3_sampmiss50_q30_snpmiss40_biallelic_mac10.bed \
    --est-bed ash_scdp3_sampmiss50_q30_snpmiss40_biallelic_mac10_thin10kbp.bed \
    --format parquet \
    --out /scratch/xe2/eh4663/LFMM/lfmm_pauci_k9_nonmf \
    --cov /scratch/xe2/eh4663/LFMM/env_cov_pauci.csv \
    --scale-cov \
    --intersect-samples \
    --threads 48
