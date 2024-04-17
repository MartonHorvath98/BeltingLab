#! bin/bash 
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
REFERENCE="data/grch38"


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

if [ -d $REFERENCE ]
then
    echo "Reference genome index already exists"
else
    # Create the reference genome index folder
    mkdir -p ${REFERENCE}

    # Create reference genom index
    source bin/makeGrch38.sh -d ${REFERENCE} -j ${CORES} -b "hisat2_index"
fi

# Salmon mapping-based quantification
echo "#########################################################"
echo "STAR alignment with salmon alignment-based quantification"
echo "#########################################################"

# Create the STAR results forlder
star_folder="results/03_star"
mkdir -p ${star_folder}
# Create the Salmon results forlder
salmon_folder="results/04_salmon"
mkdir -p ${salmon_folder}

# STAR mapping && Salmon quantification
source bin/alignSTAR.sh -i ${trim_folder} -j ${CORES} -r ${REFERENCE} -s ${star_folder} -o ${salmon_folder}

# Hisat2 alignment
echo "###############################################################"
echo "Hisat2 splice-aware alignment with featureCounts quantification"
echo "###############################################################"

# Create the Hisat2 results directory
hisat_folder="results/05_hisat"
mkdir -p ${hisat2_folder}
# Create the featureCounts results directory
count_folder="results/06_counts"
mkdir -p ${count_folder}

# Hisat2 alignment
source bin/alignHisat.sh -i ${trim_folder} -j ${CORES} -r ${ref_genome} -b "grch38_index" -o ${hisat2_folder} -c ${count_folder}

# Create final report
echo "####################################"
echo "Create final report                 "
echo "####################################"

# Create the final report directory
picard_folder="results/07_picard"
mkdir -p ${picard_folder}
report_folder="results/08_report"
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

