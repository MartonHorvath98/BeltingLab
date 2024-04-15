#! bin/bash

# Take user arguments: input-, config- and output directories
while getopts i:b:p:o: flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        b) bam_dir=${OPTARG};;
        p) picard_output=${OPTARG};;
        o) output_dir=${OPTARG};;
    esac
done

# Extract alignment metrics with picard tools from the sorted bam files
if [ -d "$bam_dir" ]
then
    # Loop through the files in the input directory
    samples=$(ls ${bam_dir}/*sorted.bam)
    while read -r sample
    do
        echo "Extracting alignment metrics with Picard on sample: $(basename $sample)"
        # Extract alignment metrics with Picard
        name=$(basename $sample)
        mkdir -p ${picard_output}/${name}
        picard CollectAlignmentSummaryMetrics I=$sample O=${picard_output}/${name}/alignment_metrics.txt
        # Extract insert size metrics with Picard
        picard CollectInsertSizeMetrics I=$sample O=${picard_output}/${name}/insert_size_metrics.txt H=${picard_output}/${name}/insert_size_histogram.pdf
        picard CollectGcBiasMetrics I=$sample O=${picard_output}/${name}/gc_bias_metrics.txt CHART=${picard_output}/${name}/gc_bias_chart.pdf S=${picard_output}/${name}/gc_summary.txt
    done <<< "$samples"
fi

# Create the multiqc report
multiqc -o $output_dir $input_dir