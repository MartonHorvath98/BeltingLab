#! bin/bash
# Marton Horvath, April 2024

# Create environment
mamba env create -f config/environment.yml

#Set up directory tree
mamba activate glioblastoma

# Run the QC pipeline
echo "############################"
echo "Subsampling raw data files #"
echo "############################"

# Input folder 
input_folder="$1"
# Remove writing permissions on the input folder
chmod a-w ${input_folder}/*

# Raw data directory
raw_data="data/00_raw"
mkdir -p $raw_data

# Create soft links to the raw data from the input folder
ln -s $input_folder/* $raw_data

# Sampled data directory
sample_data="data/01_sample"
mkdir -p $sample_data

# Subsample the raw data
source bin/getSamples.sh -i $raw_data -c "config/" -o $sample_data


# Run the QC pipeline
echo "####################################"
echo "Running quality check pipeline      "
echo "####################################"

# Create the QC results directory
qc_foder="results/01_QC"
mkdir -p $qc_folder

# Run the QC pipeline
source bin/QC.sh -d $raw_data -o $qc_folder

# Run the gene prediction pipeline
echo "####################################"
echo "Create reference Grch38 genome index"
echo "####################################"

# Create the reference genome index folder
ref_genome="grch38"
mkdir -p "${ref_genome}"
mkdir -p ${input_folder}/${ref_genome}

source bin/makeGrch38.sh -d ${input_folder}/${ref_genome} -j 4 -b "hisat2_index"
for i in $(ls ${input_folder}/${ref_genome}/*); 
do
    name=$(basename $i)
    ln -s ${input_folder}/${ref_genome}/${name} ${ref_genome}/${name}
done