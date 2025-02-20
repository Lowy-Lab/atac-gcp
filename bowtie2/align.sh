#!/bin/bash
#Manifest provided as a command line arg, and the batch task index is provided as an environment variable
#Project dir is also provided as an environment variable
dos2unix $1
sample=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" $1 | cut -d , -f 1)
path=${2}/outs/${sample}
mkdir -p ${path}/logs
#gcloud storage cp -r gs://${2}${sample}/trimmed/ .
bowtie2 -x GRCh38_noalt_decoy_as --very-sensitive -X 2000 -k 10 --no-discordant -p 8 -1 ${path}/trimmed/*_val_1.fq.gz -2 ${path}/trimmed/*_val_2.fq.gz 2> ${path}/logs/bowtie2.log | samtools view -bS -q30 - > ${path}/aligned.bam
samtools sort -@ 4 -o ${path}/aligned.bam ${path}/sorted.bam
samtools index -@ 4 ${path}/sorted.bam
gcloud storage cp ${sample}.sorted* gs://${2}${sample}/bams/
gcloud storage cp ${sample}.log gs://${2}${sample}/bams/