gcloud cp $MANIFEST samples.txt
#file=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt | cut -d , -f 2)
sample=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt | cut -d , -f 1)
gcloud storage cp -r gs://${PROJECT_DIR}/${sample}/trimmed/ .
bowtie2 -x GRCh38_noalt_decoy_as --very-sensitive -X 1000 -1 trimmed/*_1_trimmed.fq.gz -2 trimmed/*_2_trimmed.fq.gz -S ${sample}.sam \
| samtools view -bS - > ${sample}.bam
samtools sort -@ 4 -o ${sample}.sorted.bam ${sample}.bam
samtools index -@ 4 ${sample}.sorted.bam

gcloud storage cp ${sample}.sorted.b* gs://${PROJECT_DIR}/${sample}/bams/