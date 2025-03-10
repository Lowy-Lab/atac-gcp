#!/bin/bash
# Default values
input_bam=""
chrom_sizes=""
output_format="none"
threshold="4"
p_value="0.01"
q_value=""
threads=1
samplesheet=""
work_dir=""
frag_range="auto"
blacklist=""
memory="2G"

# Help function
show_help() {
    echo "Usage: $0 -i <input_bam> [-c <chrom_sizes>] [-f <output_format>] [-t <threshold>] [-p <p_value>] [-q <q_value>] [-@ <threads>] [-r <fragment_range>]"
    echo ""
    echo "Required arguments:"
    echo "  -i    Input BAM file"
    echo ""
    echo "Optional arguments:"
    echo "  -c    Chromosome sizes file"
    echo "  -f    Output format (wig, bigwig, np_array), default: wig"
    echo "  -t    Threshold value"
    echo "  -p    P-value threshold, default: 0.01"
    echo "  -q    Q-value threshold (if set, p-value is ignored)"
    echo "  -@    Number of threads to use, default: 1"
    echo "  -r    PE fragment range to call peaks on, default: auto"
    echo "  -h    Show this help message"
    echo "  -w    The working directory for gcs"
    echo "  -s    The samplesheet file path"
    echo "  -b    The encode blacklist"
    echo "  -m     The memory to use for samtools"


    exit 1
}

# Parse arguments
while getopts ":i:c:f:t:p:q:@:r:h:w:s:b:m:" opt; do
    case $opt in
        i) input_bam="$OPTARG" ;;
        c) chrom_sizes="$OPTARG" ;;
        f) output_format="$OPTARG" ;;
        t) threshold="$OPTARG" ;;
        p) p_value="$OPTARG" ;;
        q) q_value="$OPTARG" ;;
        @) threads="$OPTARG" ;;
        r) frag_range="$OPTARG" ;;
        w) work_dir="$OPTARG" ;;
        s) samplesheet="$OPTARG" ;;
        b) blacklist="$OPTARG" ;;
        m) memory="$OPTARG" ;;
        h) show_help ;;
        \?) echo "Invalid option: -$OPTARG" >&2; show_help ;;
        :) echo "Option -$OPTARG requires an argument." >&2; show_help ;;
    esac
done

# Check required arguments
if [ -z "$input_bam" ]; then
    echo "Error: Input BAM file (-i) is required."
    show_help
fi
if [ -z "$work_dir" ]; then
    echo "Error: Working directory (-w) is required."
    show_help
fi
if [ -z "$samplesheet" ]; then
    echo "Error: Samplesheet file path (-s) is required."
    show_help
fi
if [ -n "$blacklist" ]; then
    gcloud storage cp "${blacklist}" blacklist.bed
fi
gcloud storage cp "${samplesheet}" samples.txt
work_dir=${work_dir%/} #strips a trailing slash if present
dos2unix samples.txt #in case it was written in windows
sample=$(awk "NR==$((${BATCH_TASK_INDEX}+1))" samples.txt | cut -d , -f 1)
gcloud storage cp ${work_dir}/outs/per_sample_outs/${sample}/bams/${input_bam} .
samtools sort -n -@ ${threads} -m $memory -o tmp.bam ${input_bam}
if [ -z "$chrom_sizes" ]  && [ "$output_format" = "bigwig" ]; then
    echo "Chromosome sizes file (-c) is required for bigwig output. No signal output will be provided"
    output_format="none"
fi
# Validate output format
if [[ ! "$output_format" =~ ^(wig|bigwig|np_array|none)$ ]]; then
    echo "Error: Invalid output format. Must be one of: wig, bigwig, np_array, none"
    exit 1
fi
if [ -n "$chrom_sizes" ]; then
    gcloud storage cp "${chrom_sizes}" chrom.sizes
fi
# Create output directory and base filename
mkdir peaks
base_name=$(basename "$input_bam" .bam)
# Log parameters
echo "===== Peak Calling Parameters ====="
echo "Input BAM: $input_bam"
echo "Chromosome sizes: ${chrom_sizes:-Not provided}"
echo "Output format: $output_format"
echo "Threshold: ${threshold:-Default}"
echo "P-value threshold: ${p_value}"
if [ -n "$q_value" ]; then
    echo "Q-value threshold: $q_value (using this instead of p-value)"
fi
echo "Threads: $threads"
echo "Fragment range: $frag_range"
echo "Working directory: $work_dir"
echo "=================================="

# Build command parameters
cmd="fseq2 callpeak tmp.bam -cpus ${threads} -o peaks -name ${base_name} -f 0 -pe_fragment_size_range ${frag_range} -t ${threshold} -p_thr ${p_value} -pe -nfr_upper_limit 150"
if [ -n "$q_value" ]; then
    cmd+=" -q_thr ${q_value}"
fi
if [ "$output_format" != "none" ]; then
    cmd+=" -sig_format ${output_format}"
fi
if [ -n "$chrom_sizes" ]; then
    cmd+=" -chrom_size_file chrom.sizes"
fi
eval $cmd
if [ -n "$blacklist" ]; then
    bedtools subtract -A -a peaks/${base_name}_peaks.narrowPeak -b blacklist.bed > peaks/${base_name}_peaks_filtered.narrowPeak
    bedtools subtract -A -a peaks/${base_name}_summits.narrowPeak -b blacklist.bed > peaks/${base_name}_summits_filtered.narrowPeak
fi

# Calculate FRIP (Fraction of Reads in Peaks)
echo "Calculating FRIP (Fraction of Reads in Peaks)..."

# Use peaks file based on whether filtering was done
peaks_file="peaks/${base_name}_peaks.narrowPeak"
if [ -n "$blacklist" ]; then
    peaks_file="peaks/${base_name}_peaks_filtered.narrowPeak"
fi

# Generate flagstat file
samtools flagstat tmp.bam > tmp.bam.flagstat

# Count reads in peaks using intersectBed as in Nextflow pipeline
READS_IN_PEAKS=$(intersectBed -a tmp.bam -b $peaks_file -bed -c -f 0.20 | awk -F '\t' '{sum += $NF} END {print sum}')

# Extract total mapped reads from flagstat
TOTAL_MAPPED=$(grep 'mapped (' tmp.bam.flagstat | grep -v "primary" | cut -d ' ' -f 1)

# Calculate FRIP score
frip=$(awk -v reads="$READS_IN_PEAKS" -v total="$TOTAL_MAPPED" 'BEGIN {printf "%.4f", reads/total}')
echo "Total mapped reads: $TOTAL_MAPPED"
echo "Reads in peaks: $READS_IN_PEAKS"
echo "FRIP score: $frip"

# Create FRIP stats file for MultiQC
mkdir -p peaks/qc
echo -e "Sample\tReads_in_peaks\tTotal_mapped\tFRIP" > peaks/qc/${base_name}_frip_stats.txt
echo -e "${sample}\t${READS_IN_PEAKS}\t${TOTAL_MAPPED}\t${frip}" >> peaks/qc/${base_name}_frip_stats.txt

# Create FRiP file as in Nextflow
echo -e "${sample}\t${frip}" > peaks/qc/${sample}.FRiP.txt

# Create MultiQC-compatible FRIP report
echo "# plot_type: 'bargraph'" > peaks/qc/${base_name}_frip_mqc.txt
echo "# section_name: 'FRIP Scores'" >> peaks/qc/${base_name}_frip_mqc.txt
echo "# description: 'Fraction of Reads in Peaks (FRIP) score for each sample'" >> peaks/qc/${base_name}_frip_mqc.txt
echo "# pconfig:" >> peaks/qc/${base_name}_frip_mqc.txt
echo "#     id: 'frip_score_plot'" >> peaks/qc/${base_name}_frip_mqc.txt
echo "#     title: 'FRIP Scores'" >> peaks/qc/${base_name}_frip_mqc.txt
echo "#     ylab: 'FRIP Score'" >> peaks/qc/${base_name}_frip_mqc.txt
echo "#     ymax: 1" >> peaks/qc/${base_name}_frip_mqc.txt
echo "#     ymin: 0" >> peaks/qc/${base_name}_frip_mqc.txt
echo "Sample\tFRIP" >> peaks/qc/${base_name}_frip_mqc.txt
echo "${sample}\t${frip}" >> peaks/qc/${base_name}_frip_mqc.txt

# Copy output files to GCS
gcloud storage cp -r peaks/* ${work_dir}/outs/per_sample_outs/${sample}/peaks
