#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p node -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J fCounts

# Marton Horvath, APR 2024 - last modified: DEC 2024
# ------------ Start a screen environment ----- #
screen -S featureCounts
# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools samtools/1.20 subread/2.0.3

# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/

# ------------ Set global variables ------------ #
# Set the number of cores
CORES=16
eval echo "Number of cores: ${CORES}"

SAMPLES="$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')"

INPUT_DIR="${WHARF}results/03_align/"
eval echo "STAR alignments from: ${INPUT_DIR}"

REF="${WHARF}refs/ctat/"

# ----------- Check executables ---------------- #
## featureCounts
COUNTS_EXE=./featureCounts
if [ ! -x "$COUNTS_EXE" ] ; then
	if ! which featureCounts ; then
		echo "Could not find subread (featureCounts) in current directory or in PATH"
		exit 1
	else
		COUNTS_EXE=`which featureCounts`
	fi
fi

echo "###############################################"
echo "Quantification with featureCounts              "
echo "###############################################"

# Create featureCounts results directory
COUNT_DIR="${WHARF}results/04_counts/"
mkdir -p ${COUNT_DIR}

# Calculate counts with featureCounts

if [ -d "$INPUT_DIR" ] ; then
    for sample in ${SAMPLES}
    do
	if [ ! -f ${COUNT_DIR}${sample}.counts.txt ] ; then 
            # Read the trimmed forward and reverse reads
            echo "Quantifying reads with featureCounts on sample: ${sample}"
            bam=${INPUT_DIR}${sample}-Aligned.sortedByCoord.out.bam

	    $COUNTS_EXE -p --countReadPairs -s 2 -T ${CORES} -a "${REF}/ref_annot.gtf" -o "${COUNT_DIR}${sample}.counts.txt" "${bam}"
	fi
    done
fi

# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
