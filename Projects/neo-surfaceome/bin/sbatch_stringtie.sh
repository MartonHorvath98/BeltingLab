#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p node -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J Stringtie

# Marton Horvath, DEC 2024 - last modified: DEC 2024
# ------------ Start a screen environment ----- #
screen -S StringtieAssembly
# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools StringTie/2.2.1 samtools/1.20 picard/3.1.1

# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/

# ------------ Set global variables ------------ #
# Set the number of cores
CORES=16
eval echo "Number of cores: ${CORES}"

SAMPLES="$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')"
CONDITIONS="$(cat config/conditions.txt)"

INPUT_DIR="${WHARF}results/03_align/"
eval echo "STAR alignments from: ${INPUT_DIR}"

REF="${WHARF}refs/ctat/"

# ----------- Check executables ---------------- #
## 1. Stringtie 
STRINGTIE_EXE=./stringtie
if [ ! -x "$STRINGTIE_EXE" ] ; then
	if ! which stringtie ; then
		echo "Could not find stringtie in current directory or in PATH"
		exit 1
	else
		STRINGTIE_EXE=`which stringtie`
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

echo "###############################################"
echo "# Genome-guided transcriptome assembly        #"
echo "###############################################"

# Create Stringtie results directory
ISOFORM_DIR="${WHARF}results/07_isoforms/"
mkdir -p ${ISOFORM_DIR}
TMP_DIR="${ISOFORM_DIR}tmp/"
mkdir -p ${TMP_DIR}

# Transcriptome assembly with stringtie

if [ -d "$INPUT_DIR" ] ; then
    for sample in ${SAMPLES}
    do
	if [ ! -f ${ISOFORM_DIR}${sample}.gtf ] ; then 

   	    echo -e "STEP1: de novo assembly of ${sample}"
            # If it is a tumor sample without replicates, create 2 by subsampling
	    if [[ ${sample} == *"tumor"* ]] ; then
		### - copy BAM file to TMP
		if [ ! -f "${TMP_DIR}${sample}-67pc.bam" ] ; then
		    echo -e "Subsampling $sample}"
		    # Copy original BAM file
		    cp "${INPUT_DIR}${sample}-Aligned.sortedByCoord.out.bam" "${TMP_DIR}${sample}-Aligned.sortedByCoord.out.bam"
		    ${SAMTOOLS_EXE} index "${TMP_DIR}${sample}-Aligned.sortedByCoord.out.bam"
		    ### - subset 33% reads
		    java -jar $PICARD DownsampleSam \
			--INPUT "${TMP_DIR}${sample}-Aligned.sortedByCoord.out.bam" \
			--OUTPUT "${TMP_DIR}${sample}-33pc.bam" \
			--RANDOM_SEED 42 --PROBABILITY 0.33 --VALIDATION_STRINGENCY SILENT
		    ${SAMTOOLS_EXE} index "${TMP_DIR}${sample}-33pc.bam"
		    ### - subset 67% reads
		    java -jar $PICARD DownsampleSam \
			--INPUT "${TMP_DIR}${sample}-Aligned.sortedByCoord.out.bam" \
			--OUTPUT "${TMP_DIR}${sample}-67pc.bam" \
			--RANDOM_SEED 42 --PROBABILITY 0.67 --VALIDATION_STRINGENCY SILENT
		    ${SAMTOOLS_EXE} index "${TMP_DIR}${sample}-67pc.bam"
		fi

		### de novo assembly with stringtie
		${STRINGTIE_EXE} -p ${CORES} --rf  "${TMP_DIR}${sample}-Aligned.sortedByCoord.out.bam" \
		    -G "${REF}ref_annot.gtf" -o "${ISOFORM_DIR}${sample}-1.gtf"
		${STRINGTIE_EXE} -p ${CORES} --rf  "${TMP_DIR}${sample}-33pc.bam" \
		    -G "${REF}ref_annot.gtf" -o "${ISOFORM_DIR}${sample}-2.gtf"
		${STRINGTIE_EXE} -p ${CORES} --rf  "${TMP_DIR}${sample}-67pc.bam" \
		    -G "${REF}ref_annot.gtf" -o "${ISOFORM_DIR}${sample}-3.gtf"
	    else
		${STRINGTIE_EXE} -p ${CORES} --rf  "${INPUT_DIR}${sample}-Aligned.sortedByCoord.out.bam" \
		    -G "${REF}ref_annot.gtf" -o "${ISOFORM_DIR}${sample}.gtf"
	    fi
	fi	
    done

    for condition in ${CONDITIONS}
    do
	if [ ! -f "${ISOFORM_DIR}${condition}-merge.gtf" ] ; then
	    # Merge replicates
	    echo -e "STEP2: merging replicates of ${condition}"
	    ${STRINGTIE_EXE} -p ${CORES} --merge \
		"${ISOFORM_DIR}${condition}-1.gtf" \
		"${ISOFORM_DIR}${condition}-2.gtf" \
		"${ISOFORM_DIR}${condition}-3.gtf" \
		-o "${ISOFORM_DIR}${condition}-merge.gtf"
	 fi
    done

    for sample in ${SAMPLES}
    do
	if [ ! -f "${ISOFORM_DIR}${sample}-isoform.gtf" ] ; then
	    # Quantify isoforms against merged assemblies
	    echo -e "STEP3: quantifying isoforms of ${sample}"
	    
	    if [[ ${sample} == *"tumor"* ]] ; then
		condition=${sample}
	    else
		condition="$(echo -e ${sample%-*[123]})"
	    fi
	
	    ${STRINGTIE_EXE} -p ${CORES} --rf -e "${INPUT_DIR}${sample}-Aligned.sortedByCoord.out.bam" \
		-b "${ISOFORM_DIR}${sample}-ballgown.ctab" -G "${ISOFORM_DIR}${condition}-merge.gtf" \
		-o "${ISOFORM_DIR}${sample}-isoform.gtf"
	fi
    done
fi

# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
