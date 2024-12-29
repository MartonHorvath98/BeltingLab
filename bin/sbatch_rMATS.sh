#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p node -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J splice-rmats

# Marton Horvath, JUL 2024 - last modified: JUL 2024

# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools rMATS-turbo/4.1.2

# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/

 ------------ Set global variables ------------ #
# Set the number of cores
CORES=16
eval echo "Number of cores: ${CORES}"

INPUT_DATA="${WHARF}results/05_star/"
REFERENCE="${WHARF}refs/gencode/"

# rMATS results directory
AS_DIR="${WHARF}results/09_AS/"
mkdir -p ${AS_DIR}

# ------- Check command executable paths ------ #
## 1. rMATS
RMATS_EXE=./rmats.py
if [ ! -x "$RMATS_EXE" ] ; then
	if ! which rmats.py ; then
		echo "Could not find rmats.py in current directory or in PATH"
		exit 1
	else
		RMATS_EXE=`which rmats.py`
	fi
fi


# ------------ Run the algorithms ------------- #
# 
# echo -e "593-2DN,593-2DH\n593-2DN,593-3D\n593-2DN,593-tumor\n673-2DN,673-2DH\n673-2DN,673-3D\n673-2DN,673-tumor\n593-tumor,673-tumor" > "config/pairs.txt"
echo -e "593-2DN,673-2DN\n593-3D,593-tumor\n673-3D,673-tumor" > "config/pairs.txt"
 
while IFS=, read -r condition1 condition2;
do
    # Read the trimmed forward and reverse reads
    alignment1="config/${condition1}.txt"
    alignment2="config/${condition2}.txt"
    echo "Alternative splicing variants: ${condition1} - ${condition2}"
        
    # Align reads with Hisat2
    OUTPUT_DIR="${AS_DIR}${condition1}-${condition2}/"
    mkdir -p ${OUTPUT_DIR}
    $RMATS_EXE --b1 "${alignment1}" --b2 "${alignment2}" --gtf "${REFERENCE}/gencode.v37.gtf" \
    -t "paired" --libType fr-firststrand --readLength 150 \
    --bi "${REFERENCE}/star_index/" --od "${OUTPUT_DIR}" --tmp "${OUTPUT_DIR}tmp/" \
    --nthread ${CORES} --paired-stats

done < "config/pairs.txt"



