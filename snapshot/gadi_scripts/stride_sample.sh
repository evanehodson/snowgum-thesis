#!/bin/bash
#PBS -N thin_200k_stride
#PBS -l ncpus=4
#PBS -l mem=16GB
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/xe2+scratch/xe2
#PBS -P xe2

module load bcftools htslib

INPUT="/scratch/xe2/eh4663/short-read/0_pre-processing/refilter/Eregnans_centurion__shortread-bwa__longread-clair3-glnexus__merged__refiltered.vcf.gz"
OUTDIR="/scratch/xe2/eh4663/short-read/0_pre-processing/refilter_200k_subset/"
SAMPLES="/g/data/xe2/eh4663/scripts/snowgums_unix.txt"

FILENAME=$(basename "$INPUT" .vcf.gz)

# 1. Subset to sample list first, keeps all FORMAT fields
bcftools view -S "$SAMPLES" --force-samples "$INPUT" -Oz -o "${OUTDIR}${FILENAME}_samplesubset.vcf.gz"
bcftools index -t "${OUTDIR}${FILENAME}_samplesubset.vcf.gz"

# 2. Get total site count
N_SITES=$(bcftools index -n "${OUTDIR}${FILENAME}_samplesubset.vcf.gz")
STRIDE=$(( N_SITES / 200000 ))
echo "N_SITES=$N_SITES  STRIDE=$STRIDE"

# 3. Pull CHROM/POS, take every Nth site
bcftools query -f '%CHROM\t%POS\n' "${OUTDIR}${FILENAME}_samplesubset.vcf.gz" \
  | awk -v s="$STRIDE" 'NR % s == 0' > "${OUTDIR}sites_stride.tsv"

wc -l "${OUTDIR}sites_stride.tsv"   # sanity check, should be ~200k

# 4. Extract the site list, keep as plain bgzipped text
bgzip -f "${OUTDIR}sites_stride.tsv"

# 5. Extract those sites directly using bcftools view -T (works on bgzipped plain text, no index needed)
bcftools view -T "${OUTDIR}sites_stride.tsv.gz" \
  "${OUTDIR}${FILENAME}_samplesubset.vcf.gz" \
  -Oz -o "${OUTDIR}${FILENAME}_200k_stride.vcf.gz"

bcftools index -c "${OUTDIR}${FILENAME}_200k_stride.vcf.gz"