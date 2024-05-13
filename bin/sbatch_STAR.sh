#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p node -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J align-STAR

# Marton Horvath, May 2024
# ------------ Start a screen environment ----- #
screen -S STAR_alignment
# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools star/2.7.11a

# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/

# ------------ Set global variables ------------ #
# Set the number of cores
CORES=16
eval echo "Number of cores: ${CORES}"

SAMPLES="$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')"

INPUT_DIR="${WHARF}results/02_trim/"
eval echo "Trimmed samples from: ${INPUT_DIR}"

REF="${WHARF}grch38"

# ------------ Create reference genom index ---- #
if [ ! -d "${REF}" ]
then
    mkdir -p ${REF}
    # Create reference genom index
    source bin/makeGrch38.sh -d ${REF} -j ${CORES} -b "hisat2_index"
else
    echo "Reference genome already exists"
    eval echo "Reference loaction: ${REF}"
fi

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

echo "###############################################"
echo "STAR alignment                                 "
echo "###############################################"

# Create featureCounts results directory
STAR_DIR="${WHARF}results/05_star/"
mkdir -p ${STAR_DIR}

# Calculate counts with featureCounts

if [ -d "$INPUT_DIR" ]
then
    for sample in ${SAMPLES}
    do
	# Read the trimmed forward and reverse reads
    echo "Aligning reads with STAR on sample: ${sample}"
	
	fw_read=${INPUT_DIR}${sample}_R1_001_val_1.fq
	rv_read=${INPUT_DIR}${sample}_R2_001_val_2.fq

        # Align reads with STAR
        $STAR_EXE --runThreadN ${CORES}\
            --genomeDir "${REF}/star_index" \
            --readFilesIn ${fw_read} ${rv_read} \
            --readFilesCommand "gunzip -c" \
            --outReadsUnmapped None \
            --outSAMtype BAM SortedByCoordinate \
            --outBAMcompression 6 \
            --outFileNamePrefix "${STAR_DIR}${sample}-" \
            --twopassMode Basic \
            --chimOutType Junctions \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 8 \
            --chimOutJunctionFormat 1 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 100000 \
            --alignIntronMax 100000 \
            --alignSJstitchMismatchNmax 5 -1 5 5 
	
    done
fi
       
# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
