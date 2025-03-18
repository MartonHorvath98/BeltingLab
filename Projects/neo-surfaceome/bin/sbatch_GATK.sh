#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p node -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J gatk-best-practice

# Marton Horvath, OCT 2024 - last modified: DEC 2024
# ------------ Start a screen environment ----- #
screen -S variant
# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools GATK/4.3.0.0 samtools/1.20 picard/3.1.1 java/OpenJDK_17+35

# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/

 ------------ Set global variables ------------ #
# Set the number of cores
CORES=16
SUBCORES=8
eval echo "Number of cores: ${CORES}"

INPUT_DIR="${WHARF}results/05_variant/"
REFERENCE="${WHARF}refs/ctat/"

SAMPLES="$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')"

# variant results directory
GATK_DIR="${INPUT_DIR}GATK/"
mkdir -p ${GATK_DIR}

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
# ------------ Run the variant caller algorithms ------------- #
if [ -d "${INPUT_DIR}" ] ; then
    TMP_DIR="${GATK_DIR}tmp/"
    for sample in ${SAMPLES}
    do
	if [ ! -f "${INPUT_DIR}${sample}-dedup.bam" ] ; then
	    echo -e "Marking duplicates of ${sample}"
	    java -jar $PICARD MarkDuplicates \
		--INPUT "${INPUT_DIR}${sample}-SURFME.sortedByCoord.out.bam" \
		--OUTPUT "${INPUT_DIR}${sample}-dedup.bam" \
		--METRICS_FILE "${GATK_DIR}${sample}-dup-metrics.txt" \
		--ASSUME_SORT_ORDER coordinate
	fi
	# create a temporary directory
	mkdir -p "${TMP_DIR}"
	
	## GATK best practice
	if [ ! -f "${GATK_DIR}${sample}-gatk.vcf.gz" ] ; then
	    ### - add read groups 
	    name="$(echo ${sample} | tr -d '-')"
	    java -jar $PICARD AddOrReplaceReadGroups \
		--INPUT "${INPUT_DIR}${sample}-dedup.bam" \
		--OUTPUT "${TMP_DIR}${sample}-dedup-grouped.bam" \
		--SORT_ORDER coordinate --RGLB "lib" --RGPL "Illumina" \
		--RGPU "unit" --RGSM "${name}" --MAX_RECORDS_IN_RAM 10000000 \
		--COMPRESSION_LEVEL 6
	    ${SAMTOOLS_EXE} index "${TMP_DIR}${sample}-dedup-grouped.bam"
	    ### - RNA-2-DNA convention conversion 
	    gatk SplitNCigarReads --input "${TMP_DIR}${sample}-dedup-grouped.bam" \
		--max-reads-in-memory 10000000 \
		--output "${TMP_DIR}${sample}-dedup-grouped-split.bam" --reference "${REFERENCE}ref_genome.fa"
	    ${SAMTOOLS_EXE} index "${TMP_DIR}${sample}-dedup-grouped-split.bam"
	    ### - Mask known variant sites
	    gatk BaseRecalibrator \
		--input "${TMP_DIR}${sample}-dedup-grouped-split.bam" \
		--reference "${REFERENCE}ref_genome.fa" \
		--known-sites "/sw/data/GATK/hg38/dbsnp_146.hg38.vcf.gz" \
		--known-sites "/sw/data/GATK/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
		--output "${TMP_DIR}${sample}-before.table"
	    gatk ApplyBQSR \
		--reference "${REFERENCE}ref_genome.fa" \
		--input "${TMP_DIR}${sample}-dedup-grouped-split.bam" \
		--bqsr-recal-file "${TMP_DIR}${sample}-before.table" \
		--output  "${TMP_DIR}${sample}-dedup-grouped-split-bqsr.bam"
	    gatk BaseRecalibrator \
		--reference "${REFERENCE}ref_genome.fa" \
		--known-sites "/sw/data/GATK/hg38/dbsnp_146.hg38.vcf.gz" \
		--known-sites "/sw/data/GATK/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
		--output "${TMP_DIR}${sample}-after.table"
	    gatk AnalyzeCovariates \
		--before-report-file "${TMP_DIR}${sample}-before.table" \
		--after-report-file "${TMP_DIR}${sample}-after.table" \
		--plots-report-file "${GATK_DIR}${sample}-BQSR.pdf"
	    ${SAMTOOLS_EXE} index "${TMP_DIR}${sample}-dedup-grouped-split-bqsr.bam"
	    ### - Somatic SNP detection via GATK
	    gatk HaplotypeCaller \
		--input "${TMP_DIR}${sample}-dedup-grouped-split-bqsr.bam" \
		--output "${GATK_DIR}${sample}-GATK.vcf.gz" \
		--reference "${REFERENCE}ref_genome.fa" \
		--output-mode EMIT_VARIANTS_ONLY
	fi

	## Clean the TMP directory
	rm -r "${TMP_DIR}"
    done
fi






