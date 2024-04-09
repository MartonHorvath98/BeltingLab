# Integrative mRNA and Proteomic Analysis of Glioblastoma cell-surface antigens

$${\color{gray}\-\ Project\ works\ done\ for\ Mattias\ Belting's\ laboratory\ at\ Lund\ University\ \-}$$

## Introduction

> [!NOTE]
> *This section provide a brief overview of the biological background relevant to the research project and explain the biological problem being addressed, its significance in the field, and a summary of the approach or hypothesis being tested.*

Glioblastoma (GBM) is the most common and one of the most aggressive forms of primary malignant brain tumour. The current therapeutic strategy *i.e., surgical resection and radio chemotherapy* does not significantly prolong patient survival, however, low mutational burden, especially in paediatric tumour offers promising directions for developing personalized cancer treatments. Recent efforts of the group developed a platform for unbiased mapping the tumour surfaceome (TS-MAP) in glioblastoma, revealing the importance of cellular spatial organization on surfaceome diversity and identifying potential targets for antibody-drug conjugates [(Governa et al., 2022)](https://doi.org/10.1073/pnas.2114456119) 

At the same time, significant differences have been reported between the tumour surfaceome at different spatial levels within the tumour bulk. Properties of the tumour microenvironment, *e.g. hypoxic stress*, that can contribute to these dynamic changes in the expression of surfaceome and endocytosed (endocytome) proteins. The traditional analysis methodology through differential mRNA expression was found to offer low indicative prowess for actual protein availability [(Schwanhäusser et al., 2011)](https://doi.org/10.1038/nature10098), that can be improved taking a proteogenomic approach integrating tumour RNA sequencing with proteomics [(Rivero-Hinojosa et al., 2021)](https://doi.org/10.1038/s41467-021-26936-y). 

Genetic anomalies similarly play pivotal roles in tumour biology, contributing to the diversity and complexity of surface antigens. Gene fusions that can critically alter gene expression and protein function were found to be one of the major genomic abnormalities in glioblastoma. [(Shah et al., 2013)](https://doi.org/10.1186/1471-2164-14-818). By integrating these insights with the understanding of the tumour surfaceome and proteomic analyses, the project aims to uncover novel biomarkers and therapeutic targets, furthering the development of personalized treatments for glioblastoma that may include antibody-drug conjugates (ADCs) and CAR-T cells directed at cell-surface proteins.

## Project aims

By comparing mRNA expression levels across both 2D and 3D cell culture systems under normoxic and hypoxic conditions, this study aims to understand the effects of hypoxia on the expression of surfaceome and endocytome and how do the effects differ in 2D cultures and spheroids. Furthermore, the question how the mRNA expression correlates with the protein levels will be assessed by a complementary proteomic analysis ensuring a robust validation of transcriptional insights. Furthermore, the exploration of mRNA variants, including abnormal splicing patterns and unique mRNA junctions aims to identify novel biomarkers and therapeutic targets to be added to the group’s curated TS classifier (SURFME).

## Sample Metadata Summary

>[!NOTE]
> *This section provides a summary of the sample metadata used in this study including: sample source, treatments or conditions applied, and library preparation parameters... essential for interpreting the results of the bioinformatics analyses.*

### Metadata Fields

- **Sample ID**: Unique identifier for each sample.
- **Sample name**: Name of the files
- **Condition/Treatment**: Any experimental treatments or conditions applied to the sample.
- **Sequencing Type**: The type of sequencing performed (e.g., Whole Genome Sequencing, RNA-Seq).
- **Index**: The adaptors used during the library preparation.
- **Description**: Parameters of the reads - e.g.: *mean fragment length, min- and max fragment length, GC content, mean coverage, etc.*

## Setting Up the Environment

>[!NOTE]
> *To conduct the computational analyses required for this project, it's crucial to set up a consistent and reproducible environment. This section provides a guide through the process of setting up the computational environment: library tree, software dependencies and creating a virtual environment.*

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
cutadapt                  4.7             py310h4b81fae_1    bioconda
fastqc                    0.12.1               hdfd78af_0    bioconda
hisat2                    2.2.1                hdbdd923_6    bioconda
pip                       24.0               pyhd8ed1ab_0    conda-forge
python                    3.10.14         hd12c33a_0_cpython    conda-forge
salmon                    1.10.1               hecfa306_2    bioconda
samtools                  1.19.2               h50ea8bc_1    bioconda
seqtk                     1.4                  he4a0461_2    bioconda
star                      2.7.11b              h43eeafb_1    bioconda
trim-galore               0.6.10               hdfd78af_0    bioconda
```

## Raw data files

First, writing permissin for all users has been removed from the files, to avoid ever overwriting them! Then, instead of reallocating them, considering disk space and GDPR compliance, a symbolic link to the file locations was created.

```bash
chmod a-w /mnt/a/Sample*

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
    samples=$(ls -d ${input_dir}/Sample_*)
    
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

# Bioinformatical pipeline 
## I. RNA-seq from cell lines (2D) and organoids (3D) grown under normoxia and hypoxia
$${\color{gray}\-\ Compare\ gene\ expression\ under\ different\ growth\ conditions.\ \-}$$


## Quality Control
This step involves the pre-processing of the data to remove:
- adapter sequences (adapter trimming)
- low-quality reads
- uncalled bases
In this step, quality assessment is performed using the TrimGalore suite, which is a wrapper script around the popular tools FastQC and the adapter trimming algorithm Cutadapt. Cutadapt is a semi-global aligner algorithm (also called free-shift), which means that the sequences are allowed to freely shift relative to each other and differences are only penalised in the overlapping region between them. The algorithm works using unit costs (alignment score) to find the optimal overlap alignments, where positive value is assigned to matching bases and penalties are given for mismatches, inserts or deletions. 
>[!IMPORTANT]
> ***It is important to check that sequence quality is similar for all samples and discard outliers. As a general rule, read quality decreases towards the 3’ end of reads, and if it becomes too low, bases should be removed to improve mappability.  The quality and/or adapter trimming may result in very short sequences (sometimes as short as 0 bp), and since alignment programs may require sequences with a certain minimum length to avoid crashes to short fragments (in the case above, below 36 bases: --length 36) should not be considered either.***

Following a prelimnary quality assessment with FastQC, we can say that the overall quality of the bases is high, over the required treshold, however there is a high precentage of adapter contamination. Every adapter match seems to fall under the Illumina adapters' list, so the flag `--illumina` will be included in the trimming! *E.g. FastQC over-represented sequences fails for VI-3429-593-2DH-1_R1_001.fastq :*
| #Sequence                                          | Count | Percentage | Possible Source                          |
|----------------------------------------------------|-------|------------|------------------------------------------|
| GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGGAGCATCTCGTAT | 119   | 1.19       | TruSeq Adapter, Index 18 (97% over 38bp) |
| GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGGAGCATCGCGTAT | 47    | 0.47       | TruSeq Adapter, Index 18 (97% over 38bp) |
| GTCGGCATGTATTAGCTCTAGAATTACCACAGTTATCCAAGTAGGAGAGG | 16    | 0.16       | No Hit                                   |
| CTTTTACTTCCTCTAGATAGTCAAGTTCGACCGTCTTCTCAGCGCTCCGC | 14    | 0.14       | No Hit                                   |
| CCGACTTCCATGGCCACCGTCCTGCTGTCTATATCAACCAACACCTTTTC | 11    | 0.11       | No Hit                                   |

```bash
# Loop through the files in the input directory
echo "Trimming adapters with Trim Galore on sample: ${sample}"
fw_read=$(ls ${input_dir}/${sample}_R1_001.fastq) # find forward read pair
rv_read=$(ls ${input_dir}/${sample}_R2_001.fastq) # find reverse read pair
        
trim_galore ${fw_read} ${rv_read} -j 4 -q 20 --length 36 --paired --illumina -o ${output_dir} 
```
The arguments mean:
- `-j/--cores [INT]` – defines the number of  of cores to be used for trimming *[default: 1]*,
- `-q/--quality [INT]` – trims low-quality ends from reads below the threshold (as Phred score) in addition to adapter removal *[default: 20]*,
- `--length [INT]` – discards reads that become shorter than length INT because of either quality or adapter trimming,
- `--paired` – this option performs length trimming of quality/adapter/RRBS trimmed reads for paired-end files. To pass the validation test, both of a sequence pair are required to have a certain minimum length (defined by `--length`),
- `--illumina` - selects the adapter class to match by cutadapt, in this case Illumina (*first 13bp of the Illumina universal adapter 'AGATCGGAAGAGC'*)
- `--fastqc_args [ARGS]` – runs FastQC and passes down arguments in the form “arg1 arg2 etc.”, if we do not wish to pass extra arguments, FastQC in default mode can be invoked by `--fastqc` argument as well (***Either one should be called at a time!***).
- `-o/--output_dir` – is used to specify where all output will be written.

## Read mapping

Once the pre-processing and quality control steps are completed the resulting, high-quality data is ready for the read mapping or alignment step. Depending on the availability of a reference genome sequence, it can happen in one of two ways: 
1. When studying an organism with a reference genome, it is possible to infer which transcripts are expressed by mapping the reads to the reference genome (Genome mapping) or transcriptome (Transcriptome mapping). Mapping reads to the genome requires no knowledge of the set of transcribed regions or the way in which exons are spliced together. This approach allows the discovery of new, unannotated transcripts.
2. When working on an organism without a reference genome, reads need to be assembled first into longer contigs (de novo assembly). These contigs can then be considered as the expressed transcriptome to which reads are re-mapped for quantification.

>[!NOTE]
>*In this case, we follow the aprroach remains genome mapping, followed by building a more accurate, splicing-sensitive transcriptome in a genome-guided manner, and then to increase the number of mapped reads by aligning them to the assembled transcriptome. However, these steps will only happen as part of the later workflow to identify gene fusions, and mutated genes.*

There are many bioinformatics tools available to perform the alignment of short reads. Here, the one used is Hisat2, one of the most popular mapper tools, which stands for “hierarchical indexing for spliced alignment of transcripts 2”. Hisat2 provides fast and sensitive alignment, in a splice-sensitive manner, which makes it ideal for aligning RNA-samples. Besides Hisat2, read will also be aligned with quasi-mapper, Salmon. The technique employed by Salmon rapidly approximates the positions where reads might originate from the transcriptome without producing an explicit alignment of reads to the genome allowing a rapid quantification of gene expression levels.

>[!NOTE]
>*Besides the mapper algorithms, we will use Samtools for processing and analysing the data. It includes tools for file format conversion and manipulation, sorting, querying, and statistics amongst other methods.*

### Prepare reference genome

>[!IMPORTANT]
>***It is crucial to download the same reference in .fasta and .gff format from the same origin, to avoid conflicts in later steps of the analysis pipeline. Here, the GRCh38 release version v.108 (2022 Oct) of H. sapiens (human) from Ensembl is used as reference.***


