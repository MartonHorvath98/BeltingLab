#! bin/bash

# Marton Horvath, Feb 2024

# Take user arguments: input-, config- and output directories
while getopts "i:c:o:" opt; do
    case "$opt" in
        i) input_dir=$OPTARG ;;
        c) config_file=$OPTARG ;;
        o) output_dir=$OPTARG ;;
    esac
done

# Read in files from the input directory
if [ -d "$input_dir" ]
then
    # Create the output directory if it does not exist
    if [ ! -d "$output_dir" ]
    then
        mkdir -p $output_dir
    fi

    # Create the samples.tsv file if it does not exist
    if [ ! -f "$output_dir/${config_file}" ]
    then    
        echo -e "sample,read1,read2" > "$output_dir/${config_file}"
        # Loop through the files in the input directory
        samples=$(ls -d ${input_dir}Sample_*)
        while read -r sample
        do
        echo "Processing sample: $sample"
	    name=$(basename $sample)
	    name=$(sed 's/Sample_//g' <<< ${name})
        fw_read=$(ls ${sample}/*_R1_001.fastq.gz)
        rv_read=$(ls ${sample}/*_R2_001.fastq.gz)

        echo -e "$name,$fw_read,$rv_read" >> "$output_dir/${config_file}"
        done <<< "$samples"
    fi

else
    echo "Input directory does not exist"
fi
