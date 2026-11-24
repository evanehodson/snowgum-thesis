#!/bin/bash
#PBS -N beagle_pauci_single
#PBS -P xe2
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l ncpus=12
#PBS -l mem=128GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/xe2+scratch/xe2

module load java/jdk-17.0.2
cd $PBS_O_WORKDIR

CHR_NUM=${CHR_NUM:-65}

BEAGLE_JAR="/home/272/eh4663/bin/beagle.27Feb25.75f.jar"
SPLIT_DIR="/scratch/xe2/eh4663/short-read/0_pre-processing/refilter/split_chrs"
INPUT="$SPLIT_DIR/chr_${CHR_NUM}.vcf.gz"
OUTDIR="/scratch/xe2/eh4663/short-read/2_gea/2.1_beagle_impute"

echo "Processing Chromosome CM0746${CHR_NUM}.1 independently..."

java -Xmx120g -jar $BEAGLE_JAR \
    gt=$INPUT \
    excludesamples=exclude_all_except_pauciflora.txt \
    out=$OUTDIR/pauciflora_chr_${CHR_NUM}_imputed \
    nthreads=$PBS_NCPUS \
    impute=true

echo "Chromosome CM0746${CHR_NUM}.1 finished."