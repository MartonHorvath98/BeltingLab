#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p node -n 1
#SBATCH -t 6:00:00
#SBATCH -J align-hisat

# Marton Horvath, April 2024

# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools HISAT2/2.2.1 samtools/1.19 subread/2.0.3 

# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/

# ------------ Set global variables ------------ #
# Set the number of cores
CORES=16
eval echo "Number of cores: ${CORES}"

SAMPLES="$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')"

INPUT_DATA="${WHARF}results/02_trim/"
eval echo "Trimmed reads at: ${INPUT_DATA}"

# ------------ Create reference genom index ---- #
REF="${WHARF}grch38"

if [ ! -d "${REF}" ]
then
    mkdir -p ${REF}
    # Create reference genom index
    source bin/makeGrch38.sh -d ${REF} -j ${CORES} -b "hisat2_index"
else
    echo "Reference genome already exists"
    eval echo "Reference loaction: ${REF}"
fi

# ------------ Create output directories ------ #
# Hisat2 results
HISAT_DIR="${WHARF}results/03_hisat/"
mkdir -p ${HISAT_DIR}

# featureCounts results
COUNT_DIR="${WHARF}results/04_counts/"
mkdir -p ${COUNT_DIR}

# ------- Check command executable paths ------ #
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

# ------------ Run the algorithms ------------- #
# Read in files from the input directory and align reads with Hisat2
for sample in ${SAMPLES}
do
    # Read the trimmed forward and reverse reads
    echo "Aligning reads with Hisat2 on sample: ${sample}"
    fw_read=${INPUT_DATA}${sample}_R1_001_val_1.fq.gz
    rv_read=${INPUT_DATA}${sample}_R2_001_val_2.fq.gz
        
    # Align reads with Hisat2
    mkdir -p "${HISAT_DIR}${sample}"
    $HISAT2_EXE --phred33 --dta --non-deterministic -p ${CORES} --novel-splicesite-outfile "${HISAT_DIR}${sample}/novel_splicesite.txt"\
    --summary-file "${HISAT_DIR}${sample}/stats.txt" --new-summary -x "${REF}/hisat/grch38_index" -1 ${fw_read} -2 ${rv_read} |\
    $SAMTOOLS_EXE view -h -bS > "${HISAT_DIR}${sample}.bam"

    # Sort the bam file
    $SAMTOOLS_EXE sort -@ ${CORES} -o "${HISAT_DIR}${sample}.sorted.bam" "${HISAT_DIR}${sample}.bam"

    # Index the sorted bam file
    $SAMTOOLS_EXE index "${HISAT_DIR}${sample}.sorted.bam"

    # Remove the unsorted bam file
    rm "${HISAT_DIR}${sample}.bam"

done
