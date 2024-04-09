# BeltingLab
## *Project works done for Mattias Belting's laboratory at Lund University*

# Integrative mRNA and Proteomic Analysis of Glioblastoma cell-surface antigens

## Introduction

[//]: # (This section provide a brief overview of the biological background relevant to the research project and explain the biological problem being addressed, its significance in the field, and a summary of the approach or hypothesis being tested.)

Glioblastoma (GBM) is the most common and one of the most aggressive forms of primary malignant brain tumour. The current therapeutic strategy *i.e., surgical resection and radio chemotherapy* does not significantly prolong patient survival, however, low mutational burden, especially in paediatric tumour offers promising directions for developing personalized cancer treatments. Recent efforts of the group developed a platform for unbiased mapping the tumour surfaceome (TS-MAP) in glioblastoma, revealing the importance of cellular spatial organization on surfaceome diversity and identifying potential targets for antibody-drug conjugates [(Governa et al., 2022)](https://doi.org/10.1073/pnas.2114456119) 

At the same time, significant differences have been reported between the tumour surfaceome at different spatial levels within the tumour bulk. Properties of the tumour microenvironment, *e.g. hypoxic stress*, that can contribute to these dynamic changes in the expression of surfaceome and endocytosed (endocytome) proteins. The traditional analysis methodology through differential mRNA expression was found to offer low indicative prowess for actual protein availability [(Schwanhäusser et al., 2011)](https://doi.org/10.1038/nature10098), that can be improved taking a proteogenomic approach integrating tumour RNA sequencing with proteomics [(Rivero-Hinojosa et al., 2021)](https://doi.org/10.1038/s41467-021-26936-y). 

Genetic anomalies similarly play pivotal roles in tumour biology, contributing to the diversity and complexity of surface antigens. Gene fusions that can critically alter gene expression and protein function were found to be one of the major genomic abnormalities in glioblastoma. [(Shah et al., 2013)](https://doi.org/10.1186/1471-2164-14-818). By integrating these insights with the understanding of the tumour surfaceome and proteomic analyses, the project aims to uncover novel biomarkers and therapeutic targets, furthering the development of personalized treatments for glioblastoma that may include antibody-drug conjugates (ADCs) and CAR-T cells directed at cell-surface proteins.

## Project aims

By comparing mRNA expression levels across both 2D and 3D cell culture systems under normoxic and hypoxic conditions, this study aims to understand the effects of hypoxia on the expression of surfaceome and endocytome and how do the effects differ in 2D cultures and spheroids. Furthermore, the question how the mRNA expression correlates with the protein levels will be assessed by a complementary proteomic analysis ensuring a robust validation of transcriptional insights. Furthermore, the exploration of mRNA variants, including abnormal splicing patterns and unique mRNA junctions aims to identify novel biomarkers and therapeutic targets to be added to the group’s curated TS classifier (SURFME).

## Sample Metadata Summary

[//]: # (This section provides a summary of the sample metadata used in this study including: sample source, treatments or conditions applied, and library preparation parameters ... essential for interpreting the results of the bioinformatics analyses)

### Metadata Fields

- **Sample ID**: Unique identifier for each sample.
- **Sample name**: Name of the files
- **Condition/Treatment**: Any experimental treatments or conditions applied to the sample.
- **Sequencing Type**: The type of sequencing performed (e.g., Whole Genome Sequencing, RNA-Seq).
- **Index**: The adaptors used during the library preparation.
- **Description**: Parameters of the reads - e.g.: *mean fragment length, min- and max fragment length, GC content, mean coverage, etc.*

## Setting Up the Environment

[//]: # (To conduct the computational analyses required for this project, it's crucial to set up a consistent and reproducible environment. This section provides a guide through the process of setting up the computational environment: library tree, software dependencies and creating a virtual environment.)

### Directory tree
```bash
.
├── README.md
├── bin
│   ├── QC.sh
│   ├── getSamples.sh
│   └── trimReads.sh
├── config
│   ├── environment.yml
│   └── samples.csv
├── data
│   ├── 00_raw
│   ├── 01_sample
│   └── ...
├── results
│   ├── 01_QC
│   ├── 02_trim
│   └── ...
└── workflow.sh
```
### Environment prerequisits

Download and install mamba through the recommended miniforge [installation](https://github.com/conda-forge/miniforge) process.
```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```
Then, use mamba to set up the virtual environment, to run the scripts.
```bash
# Create environment
mamba env create -f config/environment.yml

#Set up directory tree
mamba activate glioblastoma
```
While the environment - called 'glioblastoma' - is created, based on the `environment.yml` file, the following programs are installed:
```bash
# packages in environment at /root/miniforge3/envs/glioblastoma:
#
# Name                    Version                   Build  Channel
cutadapt                  4.7              py39hf95cd2a_1    bioconda
fastqc                    0.12.1               hdfd78af_0    bioconda
matplotlib                3.8.3            py39hf3d152e_0    conda-forge
numpy                     1.26.4           py39h474f0d3_0    conda-forge
python                    3.9.19          h0755675_0_cpython    conda-forge
seqtk                     1.4                  he4a0461_2    bioconda
trim-galore               0.6.10               hdfd78af_0    bioconda
```

## Raw data files

First, writing permissin for all users has been removed from the files, to avoid ever overwriting them! Then, instead of reallocating them, considering disk space and GDPR compliance, a symbolic link to the file locations was created.

```bash
for file in $(cat config/samples.tsv); do 
    chmod a-w /mnt/a/${file}/*; 
done

ln -s /mnt/a/Sample* data/00_raw/
```
The raw files take up 291Gb of space, hence a subset of 10 thousand transcripts was subsampled from each using seqtk (v1.4) and the bash script `getSamples.sh`. A seed was set to make sure the same read-pairs are retained for each sample. At the same time, the names of the samples as well as the paths to the forward and reverse reads were collected in the "config/samples.csv" file to ease further analysis steps.
```bash
# Create the samples.tsv file if it does not exist
if [ ! -f "${config_dir}/samples.csv" ]
then   
    # First, add a header line
    echo -e "Sample,Forward read,Reverse read" > "${config_dir}/samples.csv"
    # Loop through the files in the input directory
    samples=$(ls -d ${input_dir}Sample_*)
    
    # Loops through the folders containing the forward and reverse reads
    while read -r sample
    do
        echo "Processing sample: $sample"
        name=$(basename $sample) # saves the sample name
        fw_read=$(ls ${sample}/*_R1_001.fastq.gz) # saves the path to the fw read
        rv_read=$(ls ${sample}/*_R2_001.fastq.gz) # saves the path to the rv read

        # Adds the name, fw read, rv read to the samples.csv file
        echo -e "$name,$fw_read,$rv_read" >> "${config_dir}/samples.csv"
    done <<< "$samples"
fi

# Pseudo-randomly subsample the reads using seqtk reading in the samples.csv file
while read -r line
do
    # Extract the sample name, forward read and reverse read
    sample=$(echo $line | cut -d ',' -f1)
    fw_read=$(echo $line | cut -d ',' -f2)
    rv_read=$(echo $line | cut -d ',' -f3)

    # subsample 10.000 reads - new files take up 3.3Mb each
    echo "Subsampling reads for sample: $sample"
    seqtk sample -s100 $fw_read 10000 > "${output_dir}/$(basename ${fw_read})"
    seqtk sample -s100 $rv_read 10000 > "${output_dir}/$(basename ${rv_read})"
done < <(tail -n +2 "${config_dir}/samples.csv") # skip the header
```


