#! bin/bash

# Take user arguments: input-, config- and output directories, mode, reference, base_name and cores
while getopts i:m:j:r:b:o: flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        m) mode=${OPTARG};;
        j) cores=${OPTARG};;
        r) reference=${OPTARG};;
        b) base_name=${OPTARG};;
        o) output_dir=${OPTARG};;
    esac
done

# Check if the mode is Hisat2
if [ $mode == "Hisat2" ]
then
    # Read in files from the input directory and align reads with Hisat2
    if [ -d "$input_dir" ]
    then
        # Create the output directory if it does not exist
        if [ ! -d "$output_dir" ]
        then
            mkdir -p $output_dir
        fi
        # Loop through the files in the input directory
        samples=$(cat "config/samples.csv" | cut -d ',' -f1 | sed 's/Sample_//g')

        while read -r sample
        do
            if [ $sample == "Sample" ]
            then
                # Skip the header line
                continue
            fi
            

            # Read the trimmed forward and reverse reads
            echo "Aligning reads with Hisat2 on sample: ${sample}"
            fw_read=$(ls ${input_dir}/${sample}_R1_001_val_1.fq)
            rv_read=$(ls ${input_dir}/${sample}_R2_001_val_2.fq)
            # Align reads with Hisat2
            mkdir -p ${output_dir}/${sample}
            hisat2 --phred33 --dta --non-deterministic -p ${cores} --novel-splicesite-outfile "${output_dir}/${sample}/novel_splicesite.txt"\
            --summary-file "${output_dir}/${sample}/stats.txt" --new-summary -x ${reference}/hisat/${base_name} -1 ${fw_read} -2 ${rv_read} |\
            samtools view -h -bS > "${output_dir}/${sample}.bam"

            # Sort the bam file
            samtools sort -@ ${cores} -o "${output_dir}/${sample}.sorted.bam" "${output_dir}/${sample}.bam"

            # Index the sorted bam file
            samtools index "${output_dir}/${sample}.sorted.bam"

            # Remove the unsorted bam file
            rm "${output_dir}/${sample}.bam"

            # Calculate readcounts with featureCounts
            featureCounts -p --countReadPairs -s 2 -f -M -O -T ${cores} -a ${reference}/Homo_sapiens.GRCh38.gff\
            -o "${output_dir}/${sample}.counts.txt" "${output_dir}/${sample}.sorted.bam"
        done <<< "${samples}"
    else
        echo "Input directory does not exist"
    fi
fi
# Check if the mode is Salmon
if [ $mode == "Salmon" ]
then
    # Read in files from the input directory and align reads with Salmon
    if [ -d "$input_dir" ]
    then
        # Create the output directory if it does not exist
        if [ ! -d "$output_dir" ]
        then
            mkdir -p $output_dir
        fi
        # Loop through the files in the input directory
        samples=$(cat "config/samples.csv" | cut -d ',' -f1 | sed 's/Sample_//g')

        while read -r sample
        do
            if [ $sample == "Sample" ]
            then
                # Skip the header line
                continue
            fi
            # Read the trimmed forward and reverse reads
            echo "Aligning reads with Salmon on sample: ${sample}"
            fw_read=$(ls ${input_dir}/${sample}_R1_001_val_1.fq)
            rv_read=$(ls ${input_dir}/${sample}_R2_001_val_2.fq)
            # Align reads with Salmon
            salmon quant -l 'ISR' -i ${reference} -1 ${fw_read} -2 ${rv_read} -p ${threads} --seqBias --validateMappings -o ${output_dir}/${sample}

            # extract libtype from the  ${output_dir}/${sample}/lib_format_counts.json file
            # LIBTYPE=$(cat ${output_dir}/${sample}/lib_format_counts.json | grep -Po '(?<="expected_format": ).*(?=,)')

        done <<< "${samples}"
    else
        echo "Input directory does not exist"
    fi
fi