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

# Check command executable paths
PICARD_EXE=./picard.jar
if [ ! -x "$PICARD_EXE" ] ; then
	if ! which picard ; then
		echo "Could not find Picard tools in current directory or in PATH"
		exit 1
	else
		PICARD_EXE=`which picard`
	fi
fi
# Extract alignment metrics with picard tools from the sorted bam files
if [ -d "$bam_dir" ]
then
    # Loop through the files in the input directory
    samples=$(ls ${bam_dir}/*sorted.bam)
    while read -r sample
    do
        echo "Extracting alignment metrics with Picard on sample: $(basename $sample)"
        # Create output folder
        name=$(basename $sample)
        mkdir -p ${picard_output}/${name}
        # Extract alignment metrics with Picard
        ${PICARD_EXE} CollectAlignmentSummaryMetrics I=$sample O=${picard_output}/${name}/alignment_metrics.txt
        # Extract insert size metrics with Picard
        ${PICARD_EXE} CollectInsertSizeMetrics I=$sample O=${picard_output}/${name}/insert_size_metrics.txt \
        H=${picard_output}/${name}/insert_size_histogram.pdf
        # Extract GC bias metrics with Picard
        ${PICARD_EXE} CollectGcBiasMetrics I=$sample O=${picard_output}/${name}/gc_bias_metrics.txt \
        CHART=${picard_output}/${name}/gc_bias_chart.pdf S=${picard_output}/${name}/gc_summary.txt
    done <<< "$samples"
fi

# Create the multiqc report
multiqc -o $output_dir $input_dir