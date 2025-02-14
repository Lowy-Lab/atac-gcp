gcloud cp $MANIFEST samples.txt
file=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt)
gcloud cp $file .
gcloud cp "${file/R1/R2}" .
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip GRCh38_noalt_as.zip