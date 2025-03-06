#!/bin/bash

# Default values
input_bam=""
logfile=""
skip_tn5=false
max_p_value=0.01
max_q_value=""
min_auc=20
min_peak_length=100
max_gap=100
peak_logfile=""
gzip_outputs=false
exclude_regions=""
work_dir=""         # New variable for working directory
samplesheet=""      # New variable for samplesheet
samtools_threads=1  # Default threads for samtools
samtools_memory="1G" # Default memory for samtools
output_name=""
verbose=false

# Function to display usage information
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -t FILE     Input BAM file (required)"
    echo "  -f FILE     Write logfile. WARNING MAY BE LARGE! (default: false)"
    echo "  -D          Skip Tn5 adjustment"
    echo "  -p FLOAT    Maximum p-value for peaks (default: 0.01)"
    echo "  -q FLOAT    Maximum q-value for peaks (ignores p-value if specified)"
    echo "  -a INT      Minimum AUC for a peak (default: 200)"
    echo "  -l INT      Minimum peak length in bp (default: 100)"
    echo "  -g INT      Maximum distance between linked significant sites (default: 100)"
    echo "  -P FILE     Call peaks from this log file instead of analyzing BAM. Path to the peak logfile from the sample directory"
    echo "  -E FILE     BED file with regions to exclude (e.g., blacklist)"
    echo "  -w DIR      Working directory for output files"
    echo "  -s FILE     Samplesheet file path"
    echo "  -z          Gzip output files"
    echo "  -@ INT      Number of threads for samtools (default: 1)"
    echo "  -m MEM      Maximum memory per thread for samtools (default: 1G)"
    echo "  -o FILE     Output file (default:    sample_genrich.narrowPeak)"
    echo "  -v          Verbose output"
    exit 1
}

# Parse command line arguments
while getopts "t:f:Dp:q:a:l:g:P:E:w:s:z@:m:o:v:" opt; do
    case $opt in
        t)
            input_bam="$OPTARG"
            ;;
        f)
            logfile="$OPTARG"
            ;;
        D)
            skip_tn5=true
            ;;
        p)
            max_p_value="$OPTARG"
            ;;
        q)
            max_q_value="$OPTARG"
            ;;
        a)
            min_auc="$OPTARG"
            ;;
        l)
            min_peak_length="$OPTARG"
            ;;
        g)
            max_gap="$OPTARG"
            ;;
        P)
            peak_logfile="$OPTARG"
            ;;
        E)
            exclude_regions="$OPTARG"
            ;;
        w)
            work_dir="$OPTARG"
            ;;
        s)
            samplesheet="$OPTARG"
            ;;
        z)
            gzip_outputs=true
            ;;
        @)
            samtools_threads="$OPTARG"
            ;;
        m)
            samtools_memory="$OPTARG"
            ;;
        o)
            output_name="$OPTARG"
            ;;
        v)
            verbose=true
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
    esac
done


# Verify required arguments
if [ -z "$input_bam" ] && [ -z "$peak_logfile" ]; then
    echo "Error: Either input BAM (-t) or peak log file (-P) is required"
    usage
fi

# Display parsed arguments (for debugging)
echo "=== Peak Calling Parameters ==="
echo "Input BAM:        $input_bam"
echo "Log file:         $logfile"
echo "Skip Tn5:         $skip_tn5"
echo "Max p-value:      $max_p_value"
echo "Max q-value:      $max_q_value"
echo "Min AUC:          $min_auc"
echo "Min peak length:  $min_peak_length"
echo "Max gap:          $max_gap"
echo "Peak log file:    $peak_logfile"
echo "Exclude regions:  $exclude_regions"
echo "Working directory: $work_dir"
echo "Samplesheet:      $samplesheet"
echo "Gzip outputs:     $gzip_outputs"
echo "==========================="

gcloud storage cp $samplesheet samples.txt
dos2unix samples.txt #in case it was written in windows
sample=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt | cut -d , -f 1)
gcloud storage cp ${work_dir}/outs/per_sample_outs/${sample}/bams/${input_bam} .
gcloud storage cp ${work_dir}/outs/per_sample_outs/${sample}/bams/${input_bam}.bai .
gcloud storage cp $exclude_regions blacklist.bed
if [ -n "$input_bam" ]; then
    echo "Position sorting bam file for Genrich!"
    samtools sort -o tmp.bam -m $samtools_memory -@ $samtools_threads -n $input_bam
    samtools index -@ $samtools_threads tmp.bam
    input_bam="tmp.bam"
fi
if [ -n "$peak_logfile" ]; then
    echo "Getting existing peak logfile!"
    gcloud storage cp ${work_dir}/outs/per_sample_outs/${sample}/logs/Genrich/${peak_logfile} .
fi
if [ -z "$output_name" ]; then
    output_name="${sample}_genrich.narrowPeak"
fi

cmd="/tools/Genrich/Genrich"

if [ -n "$input_bam" ]; then
    cmd+=" -t $input_bam"
fi

if [ -n "$logfile" ]; then
    cmd+=" -f $logfile"
fi

if $skip_tn5; then
    cmd+=" -D"
fi

if [ -n "$max_p_value" ] && [ -z "$max_q_value" ]; then
    cmd+=" -p $max_p_value"
fi

if [ -n "$max_q_value" ]; then
    cmd+=" -q $max_q_value"
fi

if [ -n "$min_auc" ]; then
    cmd+=" -a $min_auc"
fi

if [ -n "$min_peak_length" ]; then
    cmd+=" -l $min_peak_length"
fi

if [ -n "$max_gap" ]; then
    cmd+=" -g $max_gap"
fi

if [ -n "$peak_logfile" ]; then
    cmd+=" -P $peak_logfile"
fi

if [ -n "$exclude_regions" ]; then
    cmd+=" -E blacklist.bed"
fi

if $gzip_outputs; then
    cmd+=" -z"
fi
if $verbose; then
    cmd+=" -v"
fi
cmd+= " -o $output_name"

# Execute the command
echo "Executing: $cmd"
eval "$cmd"
if [ -n "$logfile" ]; then
    gcloud storage cp $logfile ${work_dir}/outs/per_sample_outs/${sample}/logs/${logfile}
fi
gcloud storage cp "$output_name" "${work_dir}/outs/per_sample_outs/${sample}/peaks/${output_name}"
