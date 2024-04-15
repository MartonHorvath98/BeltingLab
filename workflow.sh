#! bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p core -n 32
#SBATCH -t 48:00:00
#SBATCH -J rna-core

# Marton Horvath, April 2024

# ------------ Take user arguments ------------ #
# Take user arguments: input directory
if [ "$#" -ne 1 ]
then
    echo "Usage: $0 input_folder"
    exit 1
fi
# Input folder 
input_folder="$1"
# Remove writing permissions on the input folder
chmod a-w ${input_folder}/*

# ------------ Set up the environment ------------ #
# Create environment
mamba env create -f config/environment.yml
# Activate environment
mamba activate glioblastoma

# ------------ Set global variables ------------ #
# Set the number of cores
CORES=32
SAMPLES="config/samples.csv"

# ------------ Run the workflow ------------ #
# Run the QC pipeline
echo "####################################"
echo "Running quality check pipeline      "
echo "####################################"

# Raw data directory
raw_data="data/00_raw"
mkdir -p ${raw_data}

# Create soft links to the raw data from the input folder
ln -s ${input_folder}/* ${raw_data}

if [ ! -f $SAMPLES ]
then
    # Create samples.csv file
    sample_data="data/01_sample"
    mkdir -p ${sample_data}

    # Sample the raw reads
    source bin/getSamples.sh -i ${raw_data} -c ${SAMPLES} -o ${sample_data}
fi


# Create the QC results directory
qc_foder="results/01_QC"
mkdir -p $qc_folder

# Run the QC pipeline
source bin/QC.sh -d ${raw_data} -o ${qc_folder}

# Adapter trimming
echo "####################################"
echo "Adapter trimming with TrimGalore    "
echo "####################################"

# Create the adapter trimming results directory
trim_folder="results/02_trim"
mkdir -p ${trim_folder}

# Adapter trimming
source bin/trimReads.sh -i ${raw_data} -c ${SAMPLES} -p ${CORES} -o ${trim_folder}

# Create reference genom index
echo "####################################"
echo "Create reference Grch38 genome index"
echo "####################################"

# Create the reference genome index folder
ref_genome="data/grch38"
mkdir -p ${ref_genome}

# Create reference genom index
source bin/makeGrch38.sh -d ${ref_genome} -j ${CORES} -b "hisat2_index"

# Salmon mapping-based quantification
echo "####################################"
echo "Salmon mapping-based quantification "
echo "####################################"

# Create the Salmon results forlder
salmon_folder="results/03_salmon"
mkdir -p ${salmon_folder}

# Salmon mapping-based quantification
source bin/alignReads.sh -i ${trim_folder} -m "Salmon" -j ${CORES} -r ${ref_genome} -o ${salmon_folder}

# Hisat2 alignment
echo "####################################"
echo "Splice-aware alignment with Hisat2  "
echo "####################################"

# Create the Hisat2 results directory
hisat2_folder="results/04_hisat2"
mkdir -p ${hisat2_folder}

# Hisat2 alignment
source bin/alignReads.sh -i ${trim_folder} -m "Hisat2" -j ${CORES} -r ${ref_genome} -b "grch38_index" -o ${hisat2_folder}

# Create final report
echo "####################################"
echo "Create final report                 "
echo "####################################"

# Create the final report directory
picard_folder="results/05_picard"
mkdir -p ${picard_folder}
report_folder="results/06_report"
mkdir -p ${report_folder}

# Create the final report
source bin/createReport.sh -i "results/" -b ${hisat2_folder} -p ${picard_folder} -o ${report_folder}

# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"
# Deactivate the environment
mamba deactivate

