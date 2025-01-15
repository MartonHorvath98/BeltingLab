#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p node -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J snv-calling

# Marton Horvath, OCT 2024 - last modified: DEC 2024
# ------------ Start a screen environment ----- #
screen -S variant
# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools samtools/1.20 picard/3.1.1 bcftools/1.20 java/OpenJDK_17+35

# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/

# ------------ Set global variables ------------ #
# Set the number of cores
CORES=16
SUBCORES=8
eval echo "Number of cores: ${CORES}"

INPUT_DIR="${WHARF}results/05_variant/"
REFERENCE="${WHARF}refs/ctat/"

SAMPLES="$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')"

# variant results directory
OUTPUT_DIR="${INPUT_DIR}/BCFtools/"
mkdir -p ${OUTPUT_DIR}

# ------- Check command executable paths ------ #
## 1. Samtools for indexing BAM files
SAMTOOLS_EXE=./samtools
if [ ! -x "$SAMTOOLS_EXE" ] ; then
	if ! which samtools ; then
		echo "Could not find samtools in current directory or in PATH"
		exit 1
	else
		SAMTOOLS_EXE=`which samtools`
	fi
fi
## 2. BCFtools for variant calling
BCFTOOLS_EXE=./bcftools
if [ ! -x "$BCFTOOLS_EXE" ] ; then
	if ! which bcftools ; then
		echo "Could not find bcftools in current directory or in PATH"
		exit 1
	else
		BCFTOOLS_EXE=`which bcftools`
	fi
fi

# ------------ Run the variant caller algorithms ------------- #
if [ -d "${INPUT_DIR}" ] ; then
    TMP_DIR="${OUTPUT_DIR}tmp/"
    for sample in ${SAMPLES}
    do
	if [ ! -f "${INPUT_DIR}${sample}-dedup.bam" ] ; then
	    echo -e "Marking duplicates of ${sample}"
	    java -jar $PICARD MarkDuplicates \
		--INPUT "${INPUT_DIR}${sample}-SURFME.sortedByCoord.out.bam" \
		--OUTPUT "${INPUT_DIR}${sample}-dedup.bam" \
		--METRICS_FILE "${OUTPUT_DIR}${sample}-dup-metrics.txt" \
		--ASSUME_SORT_ORDER coordinate
	fi
	# create a temporary directory
	mkdir -p "${TMP_DIR}"

	## 1.) The "Old faithful": preprocessing with picard & variant calling with bcftools
	if [ ! -f "${OUTPUT_DIR}${sample}-bcftools.vcf.gz" ] ; then	    
	    ${SAMTOOLS_EXE} index "${INPUT_DIR}${sample}-dedup.bam"
	    ### - subset 25% reads
	    java -jar $PICARD DownsampleSam \
		--INPUT "${INPUT_DIR}${sample}-dedup.bam" \
		--OUTPUT "${TMP_DIR}${sample}-25pc.bam" \
		--RANDOM_SEED 42 --PROBABILITY 0.25 --VALIDATION_STRINGENCY SILENT
	    ${SAMTOOLS_EXE} index "${TMP_DIR}${sample}-25pc.bam"
	    ### - subset 50% reads
	    java -jar $PICARD DownsampleSam \
		--INPUT "${INPUT_DIR}${sample}-dedup.bam" \
		--OUTPUT "${TMP_DIR}${sample}-50pc.bam" \
		--RANDOM_SEED 42 --PROBABILITY 0.5 --VALIDATION_STRINGENCY SILENT
	    ${SAMTOOLS_EXE} index "${TMP_DIR}${sample}-50pc.bam"
	    ### - subset 75% reads
	    java -jar $PICARD DownsampleSam \
		--INPUT "${INPUT_DIR}${sample}-dedup.bam" \
		--OUTPUT "${TMP_DIR}${sample}-75pc.bam" \
		--RANDOM_SEED 42 --PROBABILITY 0.75 --VALIDATION_STRINGENCY SILENT
	    ${SAMTOOLS_EXE} index "${TMP_DIR}${sample}-75pc.bam"
	    ### - calculating coverage
	    ${BCFTOOLS_EXE} mpileup --redo-BAQ --min-BQ 30 --max-depth 10000 --per-sample-mF \
		--annotate FORMAT/DP,FORMAT/AD --threads "${CORES}" -Ou -f "${REFERENCE}ref_genome.fa" \
		"${INPUT_DIR}${sample}-dedup.bam" "${TMP_DIR}${sample}-25pc.bam" \
		"${TMP_DIR}${sample}-50pc.bam" "${TMP_DIR}${sample}-75pc.bam" |\
                ${BCFTOOLS_EXE} call --multiallelic-caller --variants-only -p 0.05 -Oz >\
                "${OUTPUT_DIR}${sample}-bcftools.vcf.gz"
	    ### - visualize possible biases
	    ${BCFTOOLS_EXE} stats --threads ${CORES} -F "${REFERENCE}ref_genome.fasta" -s - \
		"${OUTPUT_DIR}${sample}-bcftools.vcf.gz" > "${OUTPUT_DIR}${sample}.vcf.stats"
	fi

	## Clean the TMP directory
	rm -r "${TMP_DIR}"
    done
fi






