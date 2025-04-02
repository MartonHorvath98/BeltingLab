#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p node -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J splice-rmats
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

# Marton Horvath, July 2024

# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools rMATS-turbo/4.1.2
# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/

 ------------ Set global variables ------------ #
# Input directory with alignments
INPUT_DIR="${WHARF}results/04_align/"
REF="${WHARF}refs/ctat/"
INDEX="${REF}ref_genom.fa.star.idx/"

# rMATS results directory
OUTPUT_DIR="${WHARF}results/08_AS/"
mkdir -p ${OUTPUT_DIR}

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
if [ ! -f "config/bam_files.txt" ] ; then
    ls "${INPUT_DIR}"*.bam | tr "\n" "," | sed 's/.$//' > config/bam_files.txt
fi

BAMLIST="config/bam_files.txt"

srun --exclusive --ntasks=1 --cpus-per-task=16 $RMATS_EXE --b1 "${BAMLIST}" --gtf "${REF}ref_annot.gtf" \
    -t "paired" --libType fr-firststrand --readLength 150 \
    --bi "${INDEX}" --od "${OUTPUT_DIR}" --tmp "${OUTPUT_DIR}tmp/" \
    --nthread 16 --statoff --novelSS

# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit



