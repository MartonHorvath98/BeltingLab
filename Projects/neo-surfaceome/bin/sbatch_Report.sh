#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p core -n 8
#SBATCH -t 2-00:00:00
#SBATCH -J report

# Marton Horvath, Maj 2024
# ------------ Start a screen environment ----- #
screen -S report
# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools ucsc-utilities/v421 picard/3.1.1 MultiQC/1.10.1

# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/
export REF=grch38/
# ------------ Set global variables ------------ #
# Set the number of cores
CORES=8
eval echo "Number of cores: ${CORES}"

SAMPLES="$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')"

RESULTS_DIR="${WHARF}results/"
eval echo "Result directories at: ${RESULTS_DIR}"

# ------- Check command executable paths -------
## 1. UCSC
GTF_EXE=./gtfToGenePred
if [ ! -x "$GTF_EXE" ] ; then
    if ! which gtfToGenePred ; then
        echo "Could not find UCSC suite in the current directory or in PATH"
        exit 1
    else
        STAR_EXE=`which gtfToGenePred`
    fi
fi
## 2. GATK tools
MULTIQC_EXE=./multiqc
if [ ! -x "$MULTIQC_EXE" ] ; then
    if ! which multiqc ; then
	echo "Could not find MultiQC in the current direcotry or in your PATH"
	exit 1
    else
	MULTIQC_EXE=`which multiqc`
    fi
fi

# ------------ Run the workflow ------------ #
# Create final report
echo "####################################"
echo "Create final report                 "
echo "####################################"

if [ ! -f  "${WHARF}${REF}Homo_sapiens.GRCh38.gtf.refflat" ]
then
    $GTF_EXE \
	-genePredExt \
	-geneNameAsName2 \
	-ignoreGroupsWithoutExons \
	"${WHARF}${REF}Homo_sapiens.GRCh38.gtf" /dev/stdout | \
	awk 'BEGIN { OFS="\t" } { print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10 }' > "${WHARF}${REF}Homo_sapiens.GRCh38.gtf.refflat"
fi

HISAT_DIR="${WHARF}results/03_hisat/"
eval echo "HISAT2 alignment results are: ${HISAT_DIR}"
STAR_DIR="${WHARF}results/05_star/"
eval echo "STAR alignment results are: ${STAR_DIR}"

PICARD_DIR="${WHARF}results/07_picard"
mkdir -p ${PICARD_DIR}

if [ -d "$HISAT_DIR" ] && [ -d "$STAR_DIR" ] 
then 
    mkdir -p "${PICARD_DIR}/hisat"
    mkdir -p "${PICARD_DIR}/star"

    for sample in ${SAMPLES}
    do 
	if [ ! -f "${PICARD_DIR}/hisat/${sample}.RNA.metrics" ]
	then
	    # Extract alignment metrics with Picard
	    echo "Extracting Hisat alignment metrics with Picard on sample: ${sample}"

	    java -jar $PICARD CollectRnaSeqMetrics \
		--INPUT "${HISAT_DIR}${sample}.sorted.bam" \
		--OUTPUT "${PICARD_DIR}/hisat/${sample}.RNA.metrics" \
		--REF_FLAT "${WHARF}${REF}Homo_sapiens.GRCh38.gtf.refflat" \
		--STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND

	    echo "Extracting STAR alignment metrics with Picard on sample: ${sample}"

	    java -jar $PICARD CollectRnaSeqMetrics \
		--INPUT "${STAR_DIR}${sample}-Aligned.sortedByCoord.out" \
		--OUTPUT "${PICARD_DIR}/star/${sample}.RNA.metrics" \
		--REF_FLAT "${WHARF}${REF}Homo_sapiens.GRCh38.gtf.refflat" \
		--STRAND_SPECIFICITY SECOND_READ_TRANSCRIPTION_STRAND
	fi
    done
fi

REPORT_DIR="${WHARF}results/08_report"
mkdir -p ${REPORT_DIR}
multiqc -o ${REPORT_DIR} ${RESULTS_DIR}


# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
