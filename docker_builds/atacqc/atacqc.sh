#!/bin/bash
#$1 -- manifest file for samples
#$2 -- project directory
workdir=${2%/} #strips a trailing slash if present
gcloud storage cp $1 samples.txt
dos2unix samples.txt #in case it was written in windows
sample=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt | cut -d , -f 1)
gcloud storage cp ${workdir}/outs/per_sample_outs/${sample}/bams/filtered_sorted.bam .
gcloud storage cp ${workdir}/outs/per_sample_outs/${sample}/bams/filtered_sorted.bam.bai .
Rscript /fastq/atacqc.R
gcloud storage cp library_complexity.csv ${workdir}/outs/per_sample_outs/${sample}/library_complexity.csv
gcloud storage cp fragment_size_distribution.png ${workdir}/outs/per_sample_outs/${sample}/plots/