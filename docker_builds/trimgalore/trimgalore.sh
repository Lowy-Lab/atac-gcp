#!/bin/bash
#1 -- manifest file for samples
#2 -- project directory
workdir=${2%/} #strips a trailing slash if present
gcloud storage cp $1 samples.txt
echo $1
echo $2
dos2unix samples.txt
file=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt | cut -d , -f 2)
sample=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt | cut -d , -f 1)
echo $file
echo $sample
gcloud storage cp "${file}" .
gcloud storage cp "${file/R1/R2}" .
mkdir fastqc
/tools/FastQC/fastqc *.fastq.gz -o fastqc -t 6
trim_galore --paired --nextera --cores 3 --output_dir trimmed/ *R1*.fastq.gz *R2*.fastq.gz
mkdir -p trimmed/fastqc
/tools/FastQC/fastqc trimmed/*fq.gz -o trimmed/fastqc -t 6
gcloud storage cp -r fastqc ${workdir}/outs/per_sample_outs/${sample}/
gcloud storage cp -r trimmed ${workdir}/outs/per_sample_outs/${sample}/