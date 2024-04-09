#! bin/bash

# Take user arguments: input-, config- and output directories
while getopts i:c:o: flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        c) config_dir=${OPTARG};;
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
    if [ ! -f "${config_dir}/samples.csv" ]
    then    
        echo -e "Sample,Forward read,Reverse read" > "${config_dir}/samples.csv"
        # Loop through the files in the input directory
        samples=$(ls -d ${input_dir}Sample_*)
        while read -r sample
        do
        echo "Processing sample: $sample"
        name=$(basename $sample)
        fw_read=$(ls ${sample}/*_R1_001.fastq.gz)
        rv_read=$(ls ${sample}/*_R2_001.fastq.gz)

        echo -e "$name,$fw_read,$rv_read" >> "${config_dir}/samples.csv"
        done <<< "$samples"
    fi

    # Sample the reads using seqtk reading in the samples.tsv file
    while read -r line
    do
        # Extract the sample name, forward read and reverse read
        sample=$(echo $line | cut -d ',' -f1)
        fw_read=$(echo $line | cut -d ',' -f2)
        rv_read=$(echo $line | cut -d ',' -f3)

        echo "Subsampling reads for sample: $sample"
        seqtk sample -s100 $fw_read 10000 > "${output_dir}/$(basename ${fw_read})"
        seqtk sample -s100 $rv_read 10000 > "${output_dir}/$(basename ${rv_read})"
    done < <(tail -n +2 "${config_dir}/samples.csv")

else
    echo "Input directory does not exist"
fi