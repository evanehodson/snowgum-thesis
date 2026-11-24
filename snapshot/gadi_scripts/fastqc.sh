#!/bin/bash
#PBS -N fastqc_filtered
#PBS -q normal
#PBS -P xe2
#PBS -l walltime=06:00:00
#PBS -l ncpus=4
#PBS -l mem=16GB
#PBS -l storage=gdata/xe2
#PBS -l wd
#PBS -r y
#PBS -J 1-731:80%10

module load fastqc

FILELIST=/g/data/xe2/eh4663/scripts/filtered_fastq_files.txt
OUTDIR=/g/data/xe2/eh4663/fastqc_filtered_out
mkdir -p $OUTDIR

START=$PBS_ARRAY_INDEX
END=$((START + 79))

sed -n "${START},${END}p" $FILELIST | while read -r f; do
    base=$(basename "$f" .fastq.gz)
    if [ -f "${OUTDIR}/${base}_fastqc.zip" ]; then
        echo "Skipping $base, already done"
    else
        fastqc -o $OUTDIR -t 4 "$f"
    fi
done