#! bin/bash

# Take user arguments: input directory and output directory
while getopts i:o: flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        o) output_dir=${OPTARG};;
    esac
done

# Read in files from the input directory and trim illumina adapters with trim galore
if [ -d "$input_dir" ]
then
    # Create the output directory if it does not exist
    if [ ! -d "$output_dir" ]
    then
        mkdir -p $output_dir/stats
    fi

    # Loop through the files in the input directory
    samples=$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')

    while read -r sample
    do
        if [ $sample == "Sample" ]
        then
            continue
        fi

        echo "Trimming adapters with Trim Galore on sample: ${sample}"
        fw_read=$(ls ${input_dir}/${sample}_R1_001.fastq)
        rv_read=$(ls ${input_dir}/${sample}_R2_001.fastq)
        
        trim_galore ${fw_read} ${rv_read} -j 4 -q 20 --length 36 --paired --illumina -o ${output_dir} 
    done <<< "${samples}"
else
    echo "Input directory does not exist"
fi