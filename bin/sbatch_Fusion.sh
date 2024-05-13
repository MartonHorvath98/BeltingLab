#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p node -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J star-fusion

# Marton Horvath, May 2024
# ------------ Start a screen environment ----- #
screen -S Fusion
# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools star-fusion/1.10.1

# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/
export CTAT=/sw/data/CTAT_RESOURCE_LIB/2021-03/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play
# ------------ Set global variables ------------ #
SAMPLES="$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')"
CORES="16"

INPUT_DIR="${WHARF}results/05_star/"
eval echo "STAR-alignments are in: ${INPUT_DIR}"

REF="${WHARF}grch38/"

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
## 1. SALMON
FUSION_EXE=./
if [ ! -x "$FUSION_EXE" ] ; then
    if ! which STAR-Fusion ; then
        echo "Could not find STAR-Fusion in current directory or in PATH"
        exit 1
    else
        FUSION_EXE=`which STAR-Fusion`
    fi
fi

echo "###############################################"
echo "Detecting gene fusions with STAR-Fusion        "
echo "###############################################"

# Create STAR-Fusion directory
FUSION_DIR="${WHARF}results/06_fusions/"
mkdir -p ${FUSION_DIR}

# Calculate fusion events with STAR-Fusion

if [ -d "$INPUT_DIR" ]
then
    for sample in ${SAMPLES}
    do
	# Read the trimmed forward and reverse reads
        echo "Detecting gene fusions on sample: ${sample}"
	JUNCTIONS="${INPUT_DIR}${sample}-Sj.out.tab"
	
	OUTPUT_DIR="${FUSION_DIR}${sample}/"
	mkdir -p "${OUTPUT_DIR}"
	$FUSION_EXE \
	    --genome-lib-dir "${CTAT}" \
	    -J "${junction_file}" \
	    --CPU "${CORES}" \
	    --output-dir "${OUTPUT_DIR}"        
    done
fi
       
# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
