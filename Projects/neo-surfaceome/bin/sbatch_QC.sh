#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p core -n 8
#SBATCH -t 2:00:00
#SBATCH -J rna-core

# Marton Horvath, April 2024

# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools FastQC/0.11.9 

# ------------ Export paths -------------------- #
export BIN=/castor/project/home/hmarton/bin/
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/
export PATH=${BIN}:${PATH}

# ------------ Set global variables ------------ #
# Set the number of cores
CORES=8
eval echo "Number of cores: ${CORES}"

SAMPLES="$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')"

RAW_DATA="${WHARF}data/00_raw/"
eval echo "Raw data at: ${RAW_DATA}"

# ------------ Run the workflow ------------ #
# Run the QC pipeline
echo "####################################"
echo "Running quality check pipeline      "
echo "####################################"

# Create the QC results directory
QC_DIR="${WHARF}results/01_QC/"
mkdir -p ${QC_DIR}
eval echo "QC output: ${QC_DIR}"

# Run FastQC
# Loop through the files in the input directory
for sample in ${SAMPLES}
do
    #Run FastQC on files
    echo "Opening folder: ${RAW_DATA}$sample"
    for file in $(ls ${RAW_DATA}Sample_${sample}/*.fastq.gz)
    do
	if [ ! -f "${QC_DIR}${sample}_fastqc.html"]
	then
	    echo "Running FastQC on sample: $file"
            fastqc -t 6 ${file} --outdir ${QC_DIR}
	fi
    done
done

echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
