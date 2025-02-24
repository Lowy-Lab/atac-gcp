#!/bin/bash
#$1 -- manifest file for samples
#$2 -- project directory
#$3 -- number of threads for samtools
workdir=${2%/} #strips a trailing slash if present
gcloud storage cp $1 samples.txt
dos2unix samples.txt
sample=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt | cut -d , -f 1)
gcloud storage cp ${workdir}/outs/${sample}/bams/sorted.bam .
gcloud storage cp ${workdir}/outs/${sample}/bams/sorted.bam.bai .
samtools view -@ $(($3/2)) -h sorted.bam | python3 filter_mt.py - - chrM | samtools view -@ $(($3/2)) -bh -f 3 -o filtered.bam
samtools sort -@ $3 -o filtered_sorted.bam filtered.bam
samtools index -@ $3 filtered_sorted.bam
gcloud storage cp filtered_sorted.b* ${workdir}/outs/${sample}/bams/