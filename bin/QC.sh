#! bin/bash

# Take user arguments: input directory and output directory
while getopts i:o: flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) output_dir=${OPTARG};;
    esac
done

# Read in files from the input directory and run FastQC
if [ -d "$input_dir" ]
then
    
    # Loop through the files in the input directory
    samples=$(ls ${input_dir}/*fastq)
    while read -r sample
    do
        echo "Running FastQC on sample: $sample"
        fastqc -t 6 $sample -o $output_dir
    done <<< "$samples"
else
    echo "Input directory does not exist"
fi
 