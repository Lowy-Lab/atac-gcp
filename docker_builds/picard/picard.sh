#!/bin/bash
manifest=""
workdir=""
subsample=""
output_name="subsampled_filtered_sorted.bam"
threads=6
mem="3G"
declare -a exclude_samples=()
while getopts "m:d:s:o:e:r:t:" opt; do
  case $opt in
    m) manifest="$OPTARG" ;;      # manifest file
    d) workdir="${OPTARG%/}" ;;   # project directory (strip trailing slash)
    s) subsample="$OPTARG" ;;     # subsampling option
    o) output_name="$OPTARG" ;;   # output name
    e) exclude_samples+=("$OPTARG") ;; # files to exclude (can be repeated)
    r) mem="$OPTARG" ;; # mem for samtools
    t) threads="$OPTARG" ;; #threads for samtools
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done
echo "Output name: $output_name"
gcloud storage cp $manifest samples.txt
dos2unix samples.txt #in case it was written in windows
sample=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt | cut -d , -f 1)
gcloud storage cp ${workdir}/outs/per_sample_outs/${sample}/bams/filtered_sorted.bam .
gcloud storage cp ${workdir}/outs/per_sample_outs/${sample}/bams/filtered_sorted.bam.bai .
basename_sample="${output_name%.bam}"
if [ "$subsample" = "0" ]; then
    echo "Calculating subsampling frequency based on files"
    subsample=$(/tools/pandas/bin/python /fastq/get_min_complex.py $workdir $sample "${exclude_samples[@]}")
    if [ "$subsample" != "1.0" ]; then
        echo "Subsampling frequency: $subsample"
        subsample=${subsample#0.} 
        samtools view -u -h -s "42.${subsample}" filtered_sorted.bam | samtools collate -O -u - | samtools fixmate -m - - | samtools sort -@ ${threads} -m ${mem} -u - | samtools markdup -r -m s -S -f ${basename_sample}_duplicate_metrics.txt - ${output_name}
    else
        echo "No subsampling as the library for ${sample} has a greater than or equal complexity to the value!"
        samtools collate -O -u filtered_sorted.bam | samtools fixmate -m - - | samtools sort -@ ${threads} -m ${mem} -u - | samtools markdup -r -m s -S -f ${basename_sample}_duplicate_metrics.txt - ${output_name}
    fi
elif [ -n "$subsample" ]; then
    echo "Using provided subsampling value: $subsample"
    subsample=$(/tools/pandas/bin/python find_subsample.py "${workdir}/outs/per_sample_outs/${sample}/library_complexity.csv" $subsample)
    if [ "$subsample" != "1.00" ]; then
        echo "Subsampling frequency: $subsample"
        subsample=${subsample#0.} 
        samtools view  -u -h -s "42.${subsample}" filtered_sorted.bam | samtools collate -O -u - | samtools fixmate -m - - | samtools sort -@ ${threads} -m ${mem} -u - | samtools markdup -r -m s -S -f ${basename_sample}_duplicate_metrics.txt - ${output_name}
    else
        echo "No subsampling as the library has a greater than or equal complexity to the value!"
        samtools collate -O -u filtered_sorted.bam | samtools fixmate -m - - | samtools sort -@ ${threads} -m ${mem} -u - | samtools markdup -r -m s -S -f ${basename_sample}_duplicate_metrics.txt - ${output_name}
    fi
else
    echo "No subsampling!"
    samtools collate -O -u filtered_sorted.bam | samtools fixmate -m - - | samtools sort -@ ${threads} -m ${mem} -u - | samtools markdup -r -m s -S -f ${basename_sample}_duplicate_metrics.txt - ${output_name}

fi
gcloud storage cp $output_name ${workdir}/outs/per_sample_outs/${sample}/bams/
gcloud storage cp ${output_name}.bai ${workdir}/outs/per_sample_outs/${sample}/bams/
gcloud storage cp ${basename_sample}_duplicate_metrics.txt ${workdir}/outs/per_sample_outs/${sample}/logs/