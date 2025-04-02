#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p node -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J star-fusion

# Marton Horvath, May 2024
# ------------ Start a screen environment ----- #
screen -S STAR-Fusion
# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools samtools/1.20 star/2.7.8a star-fusion/1.10.1 trinity/2.14.0 jellyfish/2.3.0 python/3.9.5 pysam/0.17.0-python3.9.5

# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/
export CTAT=${WHARF}refs/ctat/

# ------------ Set global variables ------------ #
SAMPLES="$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')"
CORES=4


SEQ_DIR=${WHARF}results/02_trim/
echo -e "Trimmed reads at: ${SEQ_DIR}"

STAR_DIR=${WHARF}results/03_align/
echo -e "Junction alignments are at: ${STAR_DIR}"

# ------- Check command executable paths -------
## 1. STAR-FUSION
FUSION_EXE=./STAR-Fusion
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
FUSION_DIR="${WHARF}results/06_fusion/"
mkdir -p ${FUSION_DIR}


# Calculate fusion events with STAR-Fusion
for sample in ${SAMPLES}
do
    #Read the trimmed forward and reverse reads
    echo "Detecting gene fusions on sample: ${sample}"

    OUTPUT_DIR="${FUSION_DIR}${sample}/"
    mkdir -p "${OUTPUT_DIR}"
        
    fw_read=${SEQ_DIR}${sample}_R1_001_val_1.fq.gz
    rv_read=${SEQ_DIR}${sample}_R2_001_val_2.fq.gz
    junction=${STAR_DIR}${sample}-Chimeric.out.junction

    ${FUSION_EXE} \
	-J "${junction}" \
	--left_fq "${fw_read}" \
	--right_fq "${rv_red}" \
	--genome_lib_dir "${CTAT}" --CPU "${CORES}" \
	--FusionInspector validate --denovo_reconstruct \
	--extract_fusion_reads --examine_coding_effect \
	--output_dir "${OUTPUT_DIR}" --outTmpDir "${OUTPUT_DIR}tmp/"
done
       
# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
