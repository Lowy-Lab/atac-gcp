#!/bin/bash
#$1 -- manifest file for samples
#$2 -- project directory
#$3 -- number of threads for bowtie2
#$4 -- reference genome path
workdir=${2%/} #strips a trailing slash if present
gcloud storage cp $1 samples.txt
dos2unix samples.txt
gcloud storage cp -r $4 .
ref_name=$(basename ${4%/}) #so folder name is consistent
sample=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt | cut -d , -f 1)
gcloud storage cp ${workdir}/outs/${sample}/trimmed/*.fq.gz .
(bowtie2 -x ${ref_name}/${ref_name} --very-sensitive -X 2000 -k 10 --no-discordant --met-file bowtie2_alignment-metrics.txt -p $3 -1 *1.fq.gz -2 *2.fq.gz) 2> bowtie2.log | samtools view -bS -q30 - > aligned.bam
samtools sort -@ 6 -o sorted.bam aligned.bam
samtools index -@ 6 sorted.bam
gcloud storage cp sorted.b* ${workdir}/outs/${sample}/bams/
gcloud storage cp bowtie2.log ${workdir}/outs/${sample}/logs/bowtie2.log
gcloud storage cp bowtie2_alignment-metrics.txt ${workdir}/${sample}/logs/bowtie2_alignment-metrics.txt