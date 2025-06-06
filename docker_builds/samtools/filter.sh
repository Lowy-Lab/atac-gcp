#!/bin/bash
#$1 -- manifest file for samples
#$2 -- project directory
#$3 -- number of threads for samtools
#$4 -- memory for each thread in samtools
#Note that for the filtering step, each samtool process (view1 and view2 in this case) will use half of the threads specified in the command line
#The memory only applies to the sorting step
workdir=${2%/} #strips a trailing slash if present
gcloud storage cp $1 samples.txt
dos2unix samples.txt
sample=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt | cut -d , -f 1)
gcloud storage cp ${workdir}/outs/per_sample_outs/${sample}/bams/sorted.bam .
gcloud storage cp ${workdir}/outs/per_sample_outs/${sample}/bams/sorted.bam.bai .
chroms=($(samtools idxstats sorted.bam | head -n-1 | cut -f1 | tr '\n' ' ' | python3 get_chrom.py "human"))
samtools view -@ $(($3/2)) -h sorted.bam | python3 filter_mt.py - - "${chroms[@]}" | samtools view -@ $(($3/2)) -bh -f 3 -o filtered.bam
rm sorted.bam
rm sorted.bam.bai
samtools sort -@ $3 -m $4 -o filtered_sorted.bam filtered.bam
samtools index -@ $3 filtered_sorted.bam
gcloud storage cp filtered_sorted.b* ${workdir}/outs/per_sample_outs/${sample}/bams/