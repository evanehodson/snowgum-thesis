#!/bin/bash
#PBS -N beagle_impute
#PBS -P xe2
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l ncpus=12
#PBS -l mem=128GB
#PBS -l jobfs=10GB
#PBS -l storage=gdata/xe2+scratch/xe2

module load java/jdk-17.0.2
cd $PBS_O_WORKDIR

BEAGLE_JAR="/home/272/eh4663/bin/beagle.27Feb25.75f.jar"
INPUT="/scratch/xe2/eh4663/short-read/0_pre-processing/refilter/Eregnans_centurion__shortread-bwa__longread-clair3-glnexus__merged__refiltered.vcf.gz"
OUTDIR="/scratch/xe2/eh4663/short-read/2_gea/2.1_beagle_impute"

java -Xmx120g -jar $BEAGLE_JAR \
    gt=$INPUT \
    excludesamples=exclude_all_except_niphophila.txt \
    out=$OUTDIR/niphophila_imputed \
    nthreads=$PBS_NCPUS \
    impute=true