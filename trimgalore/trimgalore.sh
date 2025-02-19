#!/bin/bash
gcloud cp $MANIFEST samples.txt
file=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt | cut -d , -f 2)
sample=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt | cut -d , -f 1)
gcloud cp $file .
gcloud cp "${file/R1/R2}" .
trimg_galore --paired --cores 3 --output_dir trimmed/ $file "${file/R1/R2}"
gcloud storage cp -r trimmed gs://${PROJECT_DIR}/${sample}/trimmed