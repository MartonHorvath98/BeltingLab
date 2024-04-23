#! bin/bash

# Take user arguments: input-, config- and output directories, mode, reference, base_name and cores
while getopts i:j:r:b:o:c: flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        j) threads=${OPTARG};;
        r) reference=${OPTARG};;
        b) base_name=${OPTARG};;
        o) hisat_output=${OPTARG};;
        c) count_output=${OPTARG};;
    esac
done
# ------- Set global variable -------
SAMPLES=$(cat "config/samples.csv" | cut -d ',' -f1 | sed 's/Sample_//g')
# ------- Check command executable paths -------
## 1. Hisat2
HISAT2_EXE=./hisat2
if [ ! -x "$HISAT2_EXE" ] ; then
	if ! which hisat2 ; then
		echo "Could not find hisat2 in current directory or in PATH"
		exit 1
	else
		HISAT2_EXE=`which hisat2`
	fi
fi
## 2. Samtools
SAMTOOLS_EXE=./samtools
if [ ! -x "$SAMTOOLS_EXE" ] ; then
	if ! which samtools ; then
		echo "Could not find samtools in current directory or in PATH"
		exit 1
	else
		SAMTOOLS_EXE=`which samtools`
	fi
fi
## 2. featureCounts
COUNTS_EXE=./featureCounts
if [ ! -x "$COUNTS_EXE" ] ; then
	if ! which featureCounts ; then
		echo "Could not find subread (featureCounts) in current directory or in PATH"
		exit 1
	else
		COUNTS_EXE=`which featureCounts`
	fi
fi

# ------- Run the algorithms -------
# Read in files from the input directory and align reads with Hisat2
if [ -d "$input_dir" ]
then
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
        mkdir -p ${hisat_output}/${sample}
        $HISAT2_EXE --phred33 --dta --non-deterministic -p ${threads} --novel-splicesite-outfile "${hisat_output}/${sample}/novel_splicesite.txt"\
        --summary-file "${output_dir}/${sample}/stats.txt" --new-summary -x ${reference}/hisat/${base_name} -1 ${fw_read} -2 ${rv_read} |\
        $SAMTOOLS_EXE view -h -bS > "${hisat_output}/${sample}.bam"

        # Sort the bam file
        $SAMTOOLS_EXE sort -@ ${threads} -o "${hisat_output}/${sample}.sorted.bam" "${hisat_output}/${sample}.bam"

        # Index the sorted bam file
        $SAMTOOLS_EXE index "${hisat_output}/${sample}.sorted.bam"

        # Remove the unsorted bam file
        rm "${hisat_output}/${sample}.bam"

        # Calculate readcounts with featureCounts
        $COUNTS_EXE -p --countReadPairs -s 2 -T ${threads} -a ${reference}/Homo_sapiens.GRCh38.gtf\
        -o "${count_output}/${sample}.counts.txt" "${hisat_output}/${sample}.sorted.bam"
    done <<< "${SAMPLES}"
else
    echo "Input directory does not exist"
fi

