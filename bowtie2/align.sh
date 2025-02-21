#!/bin/bash
#Manifest provided as a command line arg, and the batch task index is provided as an environment variable
#Project dir is also provided as an environment variable
gcloud storage cp $1 samples.txt
mkdir ref_genome/
gcloud storage cp $4/* ref_genome/ 
dos2unix $1
sample=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt | cut -d , -f 1)
gcloud storage cp ${2}/${sample}/trimmed/*.fq.gz .
bowtie2 -x ref_genome/ --very-sensitive -X 2000 -k 10 --no-discordant -p $3 -1 *1.fq.gz -2 *2.fq.gz 2> bowtie2.log | samtools view -bS -q30 - > aligned.bam
samtools sort -@ 4 -o sorted.bam aligned.bam
samtools index -@ 4 sorted.bam
gcloud storage cp sorted.b* ${2}/${sample}/bams/
gcloud storage cp bowtie2.log ${2}/${sample}/logs/bowtie2.log