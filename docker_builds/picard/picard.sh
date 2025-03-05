#!/bin/bash
#$1 -- manifest file for samples
#$2 -- project directory
#$3 -- subsampling option
#$4 -- name of the output bams
#Subsequent arguments are files to exclude
manifest=""
workdir=""
subsample=""
output_name="subsampled_filtered_sorted.bam"
declare -a exclude_samples=()
while getopts "m:d:s:o:e:" opt; do
  case $opt in
    m) manifest="$OPTARG" ;;      # manifest file
    d) workdir="${OPTARG%/}" ;;   # project directory (strip trailing slash)
    s) subsample="$OPTARG" ;;     # subsampling option
    o) output_name="$OPTARG" ;;   # output name
    e) exclude_samples+=("$OPTARG") ;; # files to exclude (can be repeated)
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done
gcloud storage cp $manifest samples.txt
dos2unix samples.txt #in case it was written in windows
sample=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt | cut -d , -f 1)
gcloud storage cp ${workdir}/outs/${sample}/bams/filtered_sorted.bam .
gcloud storage cp ${workdir}/outs/${sample}/bams/filtered_sorted.bam.bai .
java -jar /tools/picard.jar AddOrReplaceReadGroups \
    I=filtered_sorted.bam \
    O=tmp.bam \
    RGID=flowcell1.lane1 \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=flowcell1.lane1.1 \
    RGSM=sample1
rm filtered_sorted.bam
if [ "$subsample" = "0" ]; then
    echo "Calculating subsampling frequency based on files"
    subsample=$(/tools/pandas/bin/python /tools/get_min_complex.py $workdir $sample "${exclude_samples[@]}")
    if [ "$subsample" != "1.00"]; then
        echo "Subsampling frequency: $subsample"
        subsample=${subsample#0.} 
        samtools view -h -b -s "42.${subsample}" tmp.bam > subsampled.bam
        rm tmp.bam
        samtools sort -o subsampled_sorted.bam -m 3G -@ 6 subsampled.bam
        samtools index subsampled_sorted.bam
        rm subsampled.bam
    else
        echo "No subsampling as the library for ${sample} has a greater than or equal complexity to the value!"
        mv tmp.bam subsampled_sorted.bam
    fi
elif [ -n "$subsample" ]; then
    echo "Using provided subsampling value: $subsample"
    subsample=$(/tools/pandas/bin/python find_subsample.py "${workdir}/outs/per_sample_outs/${sample}/library_complexity.csv" $subsample)
    if [ "$subsample" != "1.00"]; then
        echo "Subsampling frequency: $subsample"
        subsample=${subsample#0.} 
        samtools view -h -b -s "42.${subsample}" tmp.bam > subsampled.bam
        samtools sort -o subsampled_sorted.bam -m 3G -@ 6 subsampled.bam
        samtools index subsampled_sorted.bam
        rm subsampled.bam
    else
        echo "No subsampling as the library has a greater than or equal complexity to the value!"
        mv tmp.bam subsampled_sorted.bam
    fi
else
    echo "No subsampling!"
    mv tmp.bam subsampled_sorted.bam
fi
basename_sample="${output_name%.bam}"
java -jar /tools/picard.jar MarkDuplicates I=subsampled_sorted.bam O=deduped.bam M=${basename_sample}_picard_metrics.txt REMOVE_DUPLICATES=true
samtools sort -o $output_name -m 3G -@ 6 deduped.bam
samtools index -@ 6 $output_name
gcloud storage cp $output_name ${workdir}/outs/per_sample_outs/${sample}/bams/
gcloud storage cp ${output_name}.bai ${workdir}/outs/per_sample_outs/${sample}/bams/
gcloud storage cp ${basename_sample}_picard_metrics.txt ${workdir}/outs/per_sample_outs/${sample}/logs/




