#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p core -n 4
#SBATCH -t 6:00:00
#SBATCH -J rna-core

# Marton Horvath, APR 2024 - last modified: APR 2024

# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools TrimGalore/0.6.1

# ------------ Export paths -------------------- #
export BIN=/castor/project/home/hmarton/bin/
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/
export PATH=${BIN}:${PATH}

# ------------ Set global variables ------------ #
# Set the number of cores
CORES=4
eval echo "Number of cores: ${CORES}"

SAMPLES="$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')"

RAW_DATA="${WHARF}data/00_raw/"
eval echo "Raw data at: ${RAW_DATA}"

# ------------ Run the workflow ------------ #
## Adapter trimming
echo "####################################"
echo "Adapter trimming with TrimGalore    "
echo "####################################"

## Create the adapter trimming results directory
TRIM_DIR="${WHARF}results/02_trim/"
mkdir -p ${TRIM_DIR}
eval echo "Trimmed reads: ${TRIM_DIR}"

# Adapter trimming
for sample in ${SAMPLES}
do
    echo "Trimming adapters with Trim Galore on sample: ${sample}"
    fw_read=${RAW_DATA}Sample_${sample}/${sample}_R1_001.fastq.gz
    rv_read=${RAW_DATA}Sample_${sample}/${sample}_R2_001.fastq.gz
    
    if [ ! -f ${TRIM_DIR}${sample}_R1_001_val_1.fastq.gz ]
    then
	trim_galore ${fw_read} ${rv_read} -j 4 -q 20 --length 36 --paired --illumina --output_dir ${TRIM_DIR} --fastqc_args "--outdir ${TRIM_DIR}"
    fi
done

# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
