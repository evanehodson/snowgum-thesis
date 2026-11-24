#!/bin/bash
#PBS -N filter_lfmm2
#PBS -P xe2
#PBS -q normal
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l walltime=01:00:00
#PBS -l storage=gdata/xe2+scratch/xe2
#PBS -l wd

INFILE="/g/data/xe2/projects/paneuc2/notes/2026-07-22_pauci-gea/lfmm2_global_3spp.tsv"
OUTFILE="/scratch/xe2/eh4663/short-read/2_gea/lfmm2_sig_1e-5.tsv"

awk -F'\t' -v thresh=1e-5 '
NR==1 { print; next }
{
  sig=0
  for (i=4; i<=79; i+=3) {
    if ($i+0 < thresh) { sig=1; break }
  }
  if (sig) print
}' "$INFILE" > "$OUTFILE"