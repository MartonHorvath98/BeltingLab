#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p node -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J align-STAR

# Marton Horvath, Jun 2024 - last modified: DEC 2024
# ------------ Start a screen environment ----- #
screen -S STARalignment
# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools star/2.7.8a samtools/1.20

# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/
export CTAT=${WHARF}refs/ctat/

# ------------ Set global variables ------------ #
# Set the number of cores
CORES=16
eval echo "Number of cores: ${CORES}"

SAMPLES="$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')"
INPUT_DIR="${WHARF}results/02_trim/"
eval echo "Trimmed samples from: ${INPUT_DIR}"

REF="${CTAT}ref_genome.fa.star.idx/"
GENOME="${CTAT}ref_genome.fa"
GTF="${CTAT}ref_annot.gtf"

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
## 2. Samtools for indexing BAM files
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
echo "# Build STAR index                            #"
echo "###############################################"
# ------------ Create reference genom index ---- #
# Create STAR index directory
STAR_INDEX="${STAR_EXE} --runThreadN ${CORES} --runMode genomeGenerate --genomeDir ${REF} --genomeFastaFiles ${GENOME} --sjdbGTFfile ${GTF} --sjdbOverhang 100"
if [ -d "$REF" ] 
then
    echo "Star index to $(basename $GENOME) has been created"
else
    echo Running $STAR_INDEX
    mkdir -p ${REF}
    if $STAR_INDEX ; then
	echo "STAR index built"
    else
	echo "Index building failed; see error message"
    fi
fi

echo "###############################################"
echo "STAR alignment                                 "
echo "###############################################"
# -------- Align reads with STAR aligner
OUTPUT_DIR="${WHARF}results/03_align/"

# Align reads with STAR 2.7.8a

if [ ! -d "$OUTPUT_DIR" ] ; then
    # Create output directory
    mkdir -p ${OUTPUT_DIR}
fi
 
for sample in ${SAMPLES}
do
    # Read the trimmed forward and reverse reads
    echo "Aligning reads with STAR on sample: ${sample}"

    fw_read=${INPUT_DIR}${sample}_R1_001_val_1.fq.gz
    rv_read=${INPUT_DIR}${sample}_R2_001_val_2.fq.gz
 
    if [ ! -f "${OUTPUT_DIR}${sample}-Aligned.sortedByCoord.out.bam" ] ; then
	# Align reads with STAR
	$STAR_EXE \
	    --runThreadN "${CORES}" \
	    --twopassMode Basic \
	    --genomeDir "${REF}" \
	    --readFilesIn ${fw_read} ${rv_read} \
	    --readFilesCommand "gunzip -c" \
	    --outSAMtype BAM SortedByCoordinate \
	    --outBAMcompression 6 \
	    --outSAMstrandField intronMotif \
	    --outSAMunmapped Within \
	    --outFileNamePrefix "${OUTPUT_DIR}${sample}-" \
	    --outReadsUnmapped None \
	    --outFilterScoreMinOverLread 0 \
	    --outFilterMatchNminOverLread 0 \
	    --outFilterMatchNmin 0 \
	    --outFilterMismatchNmax 2 \
	    --alignSJstitchMismatchNmax 5 -1 5 5 \
	    --alignIntronMin 10 \
	    --alignIntronMax 100000 \
	    --alignMatesGapMax 100000 \
	    --chimOutType Junctions \
	    --chimSegmentMin 12 \
	    --chimJunctionOverhangMin 8 \
	    --chimOutJunctionFormat 1

	${SAMTOOLS_EXE} index "${OUTPUT_DIR}${sample}-Aligned.sortedByCoord.out.bam"
    fi
done

       
# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
