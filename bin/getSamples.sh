#! bin/bash

# Take user arguments: input-, config- and output directories
while getopts i:c:o: flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        c) config_file=${OPTARG};;
        o) output_dir=${OPTARG};;
    esac
done

# Read in files from the input directory
if [ -d "$input_dir" ]
then
    # Create the output directory if it does not exist
    if [ ! -d "$config_dir" ]
    then
        mkdir -p $config_dir
    fi

    # Create the samples.tsv file if it does not exist
    if [ ! -f "${config_file}" ]
    then    
        echo -e "Sample,Forward read,Reverse read" > "${config_file}"
        # Loop through the files in the input directory
        samples=$(ls -d ${input_dir}Sample_*)
        while read -r sample
        do
        echo "Processing sample: $sample"
        name=$(basename $sample)
        fw_read=$(ls ${sample}/*_R1_001.fastq.gz)
        rv_read=$(ls ${sample}/*_R2_001.fastq.gz)

        echo -e "$name,$fw_read,$rv_read" >> "${config_file}"
        done <<< "$samples"
    fi

else
    echo "Input directory does not exist"
fi