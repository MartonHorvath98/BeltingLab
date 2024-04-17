#! bin/bash

# Take user arguments: input-, config- and output directories, mode, reference, base_name and cores
while getopts i:j:r:b:o: flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        j) threads=${OPTARG};;
        r) reference=${OPTARG};;
        s) star_output=${OPTARG};;
        o) salmon_output=${OPTARG};;
    esac
done
# ------- Set global variable -------
SAMPLES=$(cat "config/samples.csv" | cut -d ',' -f1 | sed 's/Sample_//g')
# ------- Check command executable paths -------
## 1. STAR
STAR_EXE=./STAR
if [ ! -x "$STAR_EXE" ] ; then
    if ! which STAR ; then
        echo "Could not find STAR in current directory or in PATH"
        exit 1
    else
        STAR_EXE=`which STAR`
    fi
fi
## 2. Salmon
SALMON_EXE=./salmon
if [ ! -x "$SALMON_EXE" ] ; then
    if ! which salmon ; then
        echo "Could not find Salmon in current directory or in PATH"
        exit 1
    else
        SALMON_EXE=`which salmon`
    fi
fi

# ------- Run the algorithms -------
# Read in files from the input directory and align reads with Hisat2
if [ -d "$input_dir" ]
then
    # Align sequences with STAR
    while read -r sample
    do
        if [ $sample == "Sample" ]
        then
            # Skip the header line
            continue
        fi
        
        # Read the trimmed forward and reverse reads
        echo "Aligning reads with STAR on sample: ${sample}"
        fw_read=$(ls ${input_dir}/${sample}_R1_001_val_1.fq)
        rv_read=$(ls ${input_dir}/${sample}_R2_001_val_2.fq)
        # Align reads with STAR
        $STAR_EXE --runThreadN ${threads} \
            --genomeDir "${reference}/star_index" \
            --readFilesIn ${fw_read} ${rv_read} \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix "${star_output}/${sample}-" \
            --chimSegmentMin 20 \
            --chimOutType WithinBAM
    done <<< "${SAMPLES}"

    # Quantify STAR alignments with Salmon

    while read -r sample
    do
        if [ $sample == "Sample" ]
        then
            # Skip the header line
            continue
        fi
        
        # Read the trimmed forward and reverse reads
        echo "Quantifying features with Salmon on sample: ${sample}"
        alignment=$(ls ${star_output}/${sample}-Aligned.out.bam)
        # Align reads with Salmon
        $SALMON_EXE quant -l 'ISR' \
            -t "${reference}/HomHomo_sapiens.GRCh38.fa" \
            -g "${reference}/Homo_sapiens.GRCh38.gtf" \
            -a ${alignment} -p ${threads} \
            --seqBias \
            -o ${output_dir}/${sample}

        # extract libtype from the  ${output_dir}/${sample}/lib_format_counts.json file
        # LIBTYPE=$(cat ${output_dir}/${sample}/lib_format_counts.json | grep -Po '(?<="expected_format": ).*(?=,)')
    done <<< "${SAMPLES}"
else
    echo "Input directory does not exist"
fi