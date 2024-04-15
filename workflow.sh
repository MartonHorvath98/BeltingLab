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

# Adapter trimming
echo "####################################"
echo "Adapter trimming with TrimGalore    "
echo "####################################"

# Create the adapter trimming results directory
trim_folder="results/02_trim"
mkdir -p $trim_folder

# Adapter trimming
source bin/trimReads.sh -i config/samples.csv -o $trim_folder

# Create reference genom index
echo "####################################"
echo "Create reference Grch38 genome index"
echo "####################################"

# Create the reference genome index folder
ref_genome="grch38"
mkdir -p "${ref_genome}"

# Create reference genom index
source bin/makeGrch38.sh -d ${ref_genome} -j 4 -b "hisat2_index"

# mkdir -p ${input_folder}/${ref_genome}
# source bin/makeGrch38.sh -d ${input_folder}/${ref_genome} -j 4 -b "hisat2_index"
# cp -rs ${input_folder}/${ref_genome}/${name} ${ref_genome}/${name}


# Salmon mapping-based quantification
echo "####################################"
echo "Salmon mapping-based quantification "
echo "####################################"

# Create the Salmon results forlder
salmon_folder="results/03_salmon"

# Salmon mapping-based quantification
samples=$(cat "config/samples.csv" | cut -d ',' -f1 | sed 's/Sample_//g')

while read -r sample
do
    if [ $sample == "Sample" ]
    then
        # Skip the header line
        continue
    fi
       
    # Read the trimmed forward and reverse reads
    echo "Aligning reads with Hisat2 on sample: ${sample}"
    fw_read=$(ls ${input_dir}/${sample}_R1_001_val_1.fq)
    rv_read=$(ls ${input_dir}/${sample}_R2_001_val_2.fq)
    # Align reads with Hisat2
    hisat2 --phred33 --dta --non-deterministic -p 4 -x grch38/hisat/grch38_index -1 ${fw_read} -2 ${rv_read} -S ${output_dir}/${sample}.sam
done <<< "${samples}"



