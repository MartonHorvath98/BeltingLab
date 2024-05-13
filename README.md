# 1. Integrative mRNA and Proteomic Analysis of Glioblastoma cell-surface antigens

*Project works done for Mattias Belting's laboratory at Lund University*

- [1. Integrative mRNA and Proteomic Analysis of Glioblastoma cell-surface antigens](#1-integrative-mrna-and-proteomic-analysis-of-glioblastoma-cell-surface-antigens)
  - [1.1. Introduction](#11-introduction)
  - [1.2. Project aims](#12-project-aims)
  - [1.3. Sample Metadata Summary](#13-sample-metadata-summary)
    - [1.3.1. Metadata Fields](#131-metadata-fields)
  - [1.4. Setting Up the Environment](#14-setting-up-the-environment)
    - [1.4.1. Directory tree](#141-directory-tree)
    - [1.4.2. Environment prerequisits](#142-environment-prerequisits)
  - [1.4.3 Raw data files](#143-raw-data-files)
- [2. Bioinformatical pipeline](#2-bioinformatical-pipeline)
  - [2.1. Quality Control](#21-quality-control)
  - [2.2. Prepare reference genome](#22-prepare-reference-genome)
  - [2.3. Read mapping](#23-read-mapping)
    - [2.3.1 HISAT2](#231-hisat2)
    - [2.3.2 STAR](#232-star)
  - [2.4. Create report](#24-create-report)

## 1.1. Introduction

> [!NOTE]
> *This section provide a brief overview of the biological background relevant to the research project and explain the biological problem being addressed, its significance in the field, and a summary of the approach or hypothesis being tested.*

Glioblastoma (GBM) is the most common and one of the most aggressive forms of primary malignant brain tumour. The current therapeutic strategy *i.e., surgical resection and radio chemotherapy* does not significantly prolong patient survival, however, low mutational burden, especially in paediatric tumour offers promising directions for developing personalized cancer treatments. Recent efforts of the group developed a platform for unbiased mapping the tumour surfaceome (TS-MAP) in glioblastoma, revealing the importance of cellular spatial organization on surfaceome diversity and identifying potential targets for antibody-drug conjugates [(Governa et al., 2022)](https://doi.org/10.1073/pnas.2114456119) 

At the same time, significant differences have been reported between the tumour surfaceome at different spatial levels within the tumour bulk. Properties of the tumour microenvironment, *e.g. hypoxic stress*, that can contribute to these dynamic changes in the expression of surfaceome and endocytosed (endocytome) proteins. The traditional analysis methodology through differential mRNA expression was found to offer low indicative prowess for actual protein availability [(Schwanhäusser et al., 2011)](https://doi.org/10.1038/nature10098), that can be improved taking a proteogenomic approach integrating tumour RNA sequencing with proteomics [(Rivero-Hinojosa et al., 2021)](https://doi.org/10.1038/s41467-021-26936-y). 

Genetic anomalies similarly play pivotal roles in tumour biology, contributing to the diversity and complexity of surface antigens. Gene fusions that can critically alter gene expression and protein function were found to be one of the major genomic abnormalities in glioblastoma. [(Shah et al., 2013)](https://doi.org/10.1186/1471-2164-14-818). By integrating these insights with the understanding of the tumour surfaceome and proteomic analyses, the project aims to uncover novel biomarkers and therapeutic targets, furthering the development of personalized treatments for glioblastoma that may include antibody-drug conjugates (ADCs) and CAR-T cells directed at cell-surface proteins.

## 1.2. Project aims

By comparing mRNA expression levels across both 2D and 3D cell culture systems under normoxic and hypoxic conditions, this study aims to understand the effects of hypoxia on the expression of surfaceome and endocytome and how do the effects differ in 2D cultures and spheroids. Furthermore, the question how the mRNA expression correlates with the protein levels will be assessed by a complementary proteomic analysis ensuring a robust validation of transcriptional insights. Furthermore, the exploration of mRNA variants, including abnormal splicing patterns and unique mRNA junctions aims to identify novel biomarkers and therapeutic targets to be added to the group’s curated TS classifier (SURFME).

## 1.3. Sample Metadata Summary

>[!NOTE]
> *This section provides a summary of the sample metadata used in this study including: sample source, treatments or conditions applied, and library preparation parameters... essential for interpreting the results of the bioinformatics analyses.*

### 1.3.1. Metadata Fields

- **Sample ID**: Unique identifier for each sample.
- **Sample name**: Name of the files
- **Condition/Treatment**: Any experimental treatments or conditions applied to the sample.
- **Sequencing Type**: The type of sequencing performed (e.g., Whole Genome Sequencing, RNA-Seq).
- **Index**: The adaptors used during the library preparation.
- **Description**: Parameters of the reads - e.g.: *mean fragment length, min- and max fragment length, GC content, mean coverage, etc.*

## 1.4. Setting Up the Environment

>[!NOTE]
> *To conduct the computational analyses required for this project, it's crucial to set up a consistent and reproducible environment. This section provides a guide through the process of setting up the computational environment: library tree, software dependencies and creating a virtual environment.*

### 1.4.1. Directory tree
```bash
.
├── RStudio
│   ├── rnaseq
│   └── stats
├── bin
├── config
├── data
│   ├── 00_raw
│   └── 01_sample
├── grch38
│   ├── hisat
│   ├── salmon_index
│   └── star_index
└── results
    ├── 01_QC
    ├── 02_trim
    ├── 03_hisat
    ├── 04_counts
    ├── 05_star
    └── 06_fusion
```
### 1.4.2. Environment prerequisits

The sensitive nature of the patient tumor data makes it a requirement to work on the protected high-performance computing cluster of Uppmax (Uppsala Multidisciplinary Center for Advanced Computational Science), [BIANCA](https://www.uppmax.uu.se/resources/systems/the-bianca-cluster/). BIANCA is specifically designed for sensitive data providing a secure environment. It operates under strict access control, including two-factor authorization, connection only from a SUNET Internet Protocol ('IP') address, no direct web access, rather using an sftp port for up- and downloading data. These features ensure compliance with data protection regulations while handling and processing data that is subject to legal and ethical restrictions, e.g. the patient total RNA-seq and genome sequence date used in our research.

In order to access Bianca, an account is needed both for the Swedish User and Project Repository (SUPR) with the university portal (SWAMID) and for Uppmax. The user needs to join a SENS project and set up a two-factor authentication via their Uppmax account.

UPPMAX is a shared resource, and it uses a sceduling system - called SLURM - to ensure fair allocation. SLURM can be accessed either by submitting and running a batch job script via the `sbatch` command or by starting an interactive session with `interactive -A [project-code]`. The sbatch script must include Slurm directives and the commands necessary to execute the job: 
- **mandatory settings**: 
  - *account (-A)*: you must be a member to charge a project 
- **important settings**:  
  - *partition (-p)*: should the project run of complete nodes or only parts of them, on cores (default!) (1 node = 16 core(s)). The core hours scale with the number of course, which is important to keep in mind not to overshoot the 5000 core-hours allocated to a project each month.
  - *number of cores/nodes (-n)*
  - *time (-t)*: estimated runtime in core hours (hh:mm:ss). Always overestimate with ~50%, because jobs get killed when they reach the time limit, but users only get charged for the actual runtime. 

An example header:
```bash
#!/bin/bash -l
#SBATCH –A sens1234123 
#SBATCH –J test_job 
#SBATCH –p core 
#SBATCH –n 1 
#SBATCH –t 00:10:00 
 
module load some_software 
srun myapplication 
```
[!IMPORTANT]
>***When logging in, all users join the login node first. It is IMPORTANT, although not compulsory on BIANCA, to move onto a computing node (using either an interactive session or an sbatch script) to perform calculations. However, memory intensive processes will be down prioritized to not block access for others (and might never be finished).***

## 1.4.3 Raw data files

For security reasons, the hundreds of virtual project clusters on BIANCA are isolated from each other and the Internet. This makes file transfer to the server a little complicated: all files must be transfered through the wharf area of Bianca, that has access only to one project folder, named `<username>-<projid>` but nothing outside of it. The data transfer is executed using standard sftp protocol:
```bash
$ sftp -q <username>-<projid>@bianca-sftp.uppmax.uu.se:<username>-<projid>
```
The window will change to an SFTP environment, where bulk upload can be executed via the `mput`, bulk download via the `mget`commands.
```bash
sftp> mput -r <path/…/local_dir/> <path/…/remote_dir>
```
>[!IMPORTANT]
>***To avoid subtle error with the transfers, make sure you have write permissions for "owner" on the source files and directories. Furthermore, `sfpt put` cannot create directories, only use pre-existing one, so make sure that directories with matching names to the ones you want to copy are created in the remote folder prior to the transfer!*** 

The names of the samples as well as the paths to the forward and reverse reads were collected in the "config/samples.csv" file to ease further analysis steps and allow looping using the `getSamples.sh` script.

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
```

# 2. Bioinformatical pipeline 
*Compare gene expression RNA-seq from cell lines (2D) and organoids (3D) grown under normoxia and hypoxia*

## 2.1. Quality Control
This step involves the pre-processing of the data to remove:

- adapter sequences (adapter trimming)
- low-quality reads
- uncalled bases

In this step, quality assessment is performed using the TrimGalore suite (v0.6.1), which is a wrapper script around the popular tools FastQC and the adapter trimming algorithm Cutadapt. Cutadapt is a semi-global aligner algorithm (also called free-shift), which means that the sequences are allowed to freely shift relative to each other and differences are only penalised in the overlapping region between them. The algorithm works using unit costs (alignment score) to find the optimal overlap alignments, where positive value is assigned to matching bases and penalties are given for mismatches, inserts or deletions.

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
fw_read=${RAW_DATA}Sample_${sample}/${sample}_R1_001.fastq.gz
rv_read=${RAW_DATA}Sample_${sample}/${sample}_R2_001.fastq.gz
    
trim_galore ${fw_read} ${rv_read} -j 4 -q 20 --length 36 --paired --illumina --output_dir ${TRIM_DIR} --fastqc_args "--outdir ${TRIM_DIR}"
```
The arguments mean:
- `-j/--cores [INT]` – defines the number of  of cores to be used for trimming *[default: 1]*,
- `-q/--quality [INT]` – trims low-quality ends from reads below the threshold (as Phred score) in addition to adapter removal *[default: 20]*,
- `--length [INT]` – discards reads that become shorter than length INT because of either quality or adapter trimming,
- `--paired` – this option performs length trimming of quality/adapter/RRBS trimmed reads for paired-end files. To pass the validation test, both of a sequence pair are required to have a certain minimum length (defined by `--length`),
- `--illumina` - selects the adapter class to match by cutadapt, in this case Illumina (*first 13bp of the Illumina universal adapter 'AGATCGGAAGAGC'*)
- `--fastqc_args [ARGS]` – runs FastQC and passes down arguments in the form “arg1 arg2 etc.”, if we do not wish to pass extra arguments, FastQC in default mode can be invoked by `--fastqc` argument as well (***Either one should be called at a time!***).
- `-o/--output_dir` – is used to specify where all output will be written.

## 2.2. Prepare reference genome

Once the pre-processing and quality control steps are completed the resulting, high-quality data is ready for the read mapping or alignment step. Depending on the availability of a reference genome sequence, it can happen in one of two ways: 
1. When studying an organism with a reference genome, it is possible to infer which transcripts are expressed by mapping the reads to the reference genome (Genome mapping) or transcriptome (Transcriptome mapping). Mapping reads to the genome requires no knowledge of the set of transcribed regions or the way in which exons are spliced together. This approach allows the discovery of new, unannotated transcripts.
2. When working on an organism without a reference genome, reads need to be assembled first into longer contigs (de novo assembly). These contigs can then be considered as the expressed transcriptome to which reads are re-mapped for quantification.

>[!NOTE]
>*In this case, we follow the aprroach remains genome mapping, followed by building a more accurate, splicing-sensitive transcriptome in a genome-guided manner, and then to increase the number of mapped reads by aligning them to the assembled transcriptome. However, these steps will only happen as part of the later workflow to identify gene fusions, and mutated genes.*

>[!IMPORTANT]
>***It is crucial to download the same reference in .fasta and .gff format from the same origin, to avoid conflicts in later steps of the analysis pipeline. Here, the GRCh38 release version v.108 (2022 Oct) of H. sapiens (human) from Ensembl is used as reference.***

There are many bioinformatics tools available to perform the alignment of short reads. Here, the ones used are **Hisat2 (v.2.2.1)** (*“hierarchical indexing for spliced alignment of transcripts 2”*) and **STAR (v.2.7.11a)** (*Spliced Transcripts Alignment to a Reference*). Both Hisat2 and STAR provide sensitive and highliy efficient alignment for non-contiguous transcripts (splice-sensitive), which makes them ideal for aligning RNA-samples. In the later part of the analysis, results from the Hisat2 alignment will be used for differential expression analysis and STAR will be particularly useful for the detection of splice junctions and different isoforms.

Indexes for the Hisat2 and STAR aligners are created using the `makeGrch38.sh` script. The hisat2-build command will create eight files with the extensions *.1.ht2, .2.ht2, .3.ht2, .4.ht2, .5.ht2, .6.ht2, .7.ht2* and *.8.ht2*. These files together are the index for hisat2. Because the plans include heavy focus on alternative splicing and gene fusion events, the genome index is created with transcript annotations using the `--ss` and `--exon` flags. These take into account the prepared exons and splice sites, in HISAT2's own format (three or four tab-separated columns). The STAR index of the reference genome compirses of the binary genome sequence, suffix arrays, junction coordinates, and transcript/gene information. The suffix array method used for the construction enables rapid and memory-efficient searches, **BUT!** the file system needs to have at least 100 Gb of disk space available for the human genome (*final size ~ 30Gb*). The genome indeces for Hisat2 and STAR are generated with the following steps:
```bash
# ------- Download reference from ENSEMBL -------
# The Ensembl release version and base URLs for downloads
ENSEMBL_RELEASE=108
ENSEMBL_GENOME=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna
ENSEMBL_GFF3_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gff3/homo_sapiens/
ENSEMBL_TRANSCRIPTOME=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/cds

# ------- Create Hisat2 index -------
# Convert GFF file to GTF 
gffread <path/…/ref.gff> -T -o <path/…/ref.gtf>
# Extract splice sites and exons
hisat2_extract_splice_sites.py -v <path/…/ref.gtf> > <path/…/splicesite.tsv>
hisat2_extract_exons.py -v <path/…/ref.gtf> > <path/…/exons.tsv>
# Build HISAT2 index with transcripts annotations
hisat-build -p ${threads} --seed 100 <path/…/ref.fasta> --ss <path/…/splicesite.tsv> --exon  <path/…/exons.tsv>  <path/…/HISAT/base_name>

# ------- Create STAR index -------
# Building STAR index
star --runThreadN ${threads} --runMode genomeGenerate --genomeDir <path/…/STAR_index/> --genomeFastaFiles <path/…/ref.fa>  --sjdbGTFfile <path/…/ref.gtf> --sjdbOverhang 100
```

## 2.3. Read mapping

>[!NOTE]
>*Besides the mapper algorithms, Samtools is used for processing and analysing the data. It includes tools for file format conversion and manipulation, sorting, querying, and statistics amongst other methods.*

### 2.3.1 HISAT2

We align the RNA-seq reads to a genome in order to identify exon-exon splice junctions with the help of the `sbatch_Hisat.sh` script. 
```bash
# Read the trimmed forward and reverse reads
    echo "Aligning reads with Hisat2 on sample: ${sample}"
    fw_read=${INPUT_DATA}${sample}_R1_001_val_1.fq.gz
    rv_read=${INPUT_DATA}${sample}_R2_001_val_2.fq.gz
        
    # Align reads with Hisat2
    mkdir -p "${HISAT_DIR}${sample}"
    $HISAT2_EXE --phred33 --dta --non-deterministic -p ${CORES} --novel-splicesite-outfile "${HISAT_DIR}${sample}/novel_splicesite.txt" --summary-file "${HISAT_DIR}${sample}/stats.txt" --new-summary -x "${REF}/hisat/grch38_index" -1 ${fw_read} -2 ${rv_read} |\
    $SAMTOOLS_EXE view -h -bS > "${HISAT_DIR}${sample}.bam"

    # Sort the bam file
    $SAMTOOLS_EXE sort -@ ${CORES} -o "${HISAT_DIR}${sample}.sorted.bam" "${HISAT_DIR}${sample}.bam"

    # Index the sorted bam file
    $SAMTOOLS_EXE index "${HISAT_DIR}${sample}.sorted.bam"

    # Remove the unsorted bam file
    rm "${HISAT_DIR}${sample}.bam"
```
During the alignment step the applied parameters mean:
- `--phred33` – defines the input quality (*default: Phred+33*),
- `--dta` – (=downstream-transcriptome-assembly) report alignments tailored for transcript assemblers including StringTie,
- `--p [INT]` – sets INT number of cores to be used during the calculations,
- `--non-deterministic` – Hisat2 re-initializes its pseudo-random generator for each read using the current time, meaning that Hisat2 will not necessarily report the same alignment for two identical reads. This might be more appropriate in situations where the input consists of many identical reads,
- `--rna-strandedness 'RF'` – is used to define the strand-specificity of the RNA-seq protocol,
- `--x` – the path to the index for the reference genome including the basename of the files excluding the *.1.ht2 / etc.* extension,
- `--1` and `--2` – are comma-separated list of files containing mate 1s and mate 2s,
- `--summary-file` – sets the path to the summary file,
- `--new-summary` – makes the summary file readable by MultiQC for downstream analysis.

As, following quantification, the next downstream step will most likely be a genome-guided transcriptome assembly, the `--novel-splicesite-outfile` parameter is also set to reports a list of splice sites in the file: chromosome name tab genomic position of the flanking base on the left side of an intron tab genomic position of the flanking base on the right tab strand (+, -, and .) ‘.’ indicates an unknown strand for non-canonical splice sites. 

Following the alingment step, **Samtools (v.1.19)**, channeled through the piping operator '|', instantly converts the .sam files produced by the aligner to their binary counterpart .bam, which limits space requirements as the output file size shrinks considerably via `samtools view -h -bS > <path/…/file.bam>`. Then using the `samtools sort` and `index` command the bam files are sorted along the reference, and BAI-format index is generated.

Finally, we quantify the successfully aligned reads via **featureCounts (v.2.0.3)**, using the `sbatch_fCounts.sh` script. FeatureCounts, that is part of the subread software package, takes the resulting BAM files on its input and the genome annotation file in GFF format containing the chromosomal coordinates of features. It outputs numbers of reads assigned to features or meta-features. Each entry in the GFF annotation file is considered a feature (e.g. an exon), while a meta-feature is the aggregation of a set of features (e.g. a gene). When summarising reads at meta-feature level, read counts obtained for features included in the same meta-feature will be added up to yield the read count for the corresponding meta-feature. It also outputs stat info for the overall summarization results, including number of successfully assigned reads and number of reads that failed to be assigned due to various reasons (these reasons are included in the stat info). 
```bash
# Read the trimmed forward and reverse reads
  echo "Quantifying reads with featureCounts on sample: ${sample}"
  bam=${INPUT_DIR}${sample}.sorted.bam

	$COUNTS_EXE -p --countReadPairs -s 2 -T ${CORES} -a "${REF}/Homo_sapiens.GRCh38.gtf" -o "${COUNT_DIR}${sample}.counts.txt" "${bam}"
    
```
While counting how many reads align to each gene in a genome annotation the set parameters are:
- `-p` – specifies that input data contain paired-end reads and performs fragment counting (ie. counting read pairs), the `--countReadPairs` parameter should also be specified in addition to this parameter
- `-s [INT]` – determines the strand-specificity of the read counting (*here: 2 (reversely stranded)*),
- `-T [INT]` – specifies the number of threads to be used,
- `-a` – is a mandatory argument, the path to the annotation file,
- `-o` –  is the path to and name of the output file that includes the read counts and a separate file with summary statistics.

### 2.3.2 STAR

We align the RNA-seq reads to a genome with STAR in order to identify abnormal splicing events and chimeras (fusion genes) with the help of the `sbatch_STAR.sh` script. Following the alignment we will use **STAR-fusion (v1.10.1)** to identify candidate fusion transcripts supported by the Illumina reads. 
```bash
# Read the trimmed forward and reverse reads
echo "Aligning reads with STAR on sample: ${sample}"
	
fw_read=${INPUT_DIR}${sample}_R1_001_val_1.fq.gz
rv_read=${INPUT_DIR}${sample}_R2_001_val_2.fq.gz

# Align reads with STAR
STAR --runThreadN ${CORES} \ 
  --genomeDir "${REF}/star_index" \ 
  --readFilesIn ${fw_read} ${rv_read} \ 
  --readFilesCommand "gunzip -c" \ 
  --outReadsUnmapped None \ 
  --outSAMtype BAM SortedByCoordinate \ 
  --outBAMcompression 6 \ 
  --outFileNamePrefix "${STAR_DIR}${sample}-" \ 
  --twopassMode Basic \ 
  --chimOutType Junctions \  # **essential** to create the Chimeric.junction.out file for STAR-Fusion
  --chimOutJunctionFormat 1 \  # **essential** includes required metadata in Chimeric.out.junction file.
  --chimSegmentMin 12 \  # **essential to invoke chimeric read detection & reporting**
  --chimJunctionOverhangMin 8 \ 
  --alignSJDBoverhangMin 10 \ 
  --alignMatesGapMax 100000 \ 
  --alignIntronMax 100000 \ 
  --alignSJstitchMismatchNmax 5 -1 5 5
```
During the alignment step the applied parameters mean:
- `--runThreadN [INT]` – number of cores,
- `--genomeDir` – path to the reference STAR index,
- `--readFilesIn [R1] [R2]` – mate 1s and mate 2s reads,
- `--readFilesCommand ["gunzip -c"]` – STAR cannot work with compressed input files, so we have to create temporary unzipped versions of the trimmed fastq files
- `--outSAMtype 'BAM SortedByCoordinate'` – makes sure that the output is in sorted BAM format,
- `--outBAMcompression [INT]` – defines the compression rate of the output files (on a 0 to 10 scale, *default:6*),
- `--outFileNamePrefix` – prefix to indetify unique outputs,
- `--twopassMode Basic` – enambles the most sensitive novel splice junction discovery: during the 1st pass STAR maps with the usual parameters, then in hte 2nd pass it collects the previously detected junctions and uses the mas "annotated" junctions, what allows to detect more spliced reads mapping to novel junctions, 
- `--chimOutType Junctions` – outputs the chimeric alignments to a separate file, called *"Chimeric.out.junction"*, that is the input file for STAR-Fusion,
- `--chimOutJunctionFormat 1` – is essential to include required metadata in Chimeric.junction.out file for STAR-Fusion,
- `--chimSegmentMin [INT]` – detect fusions that map to two chromosomes, with a minimum length above [INT] on either chromosome, this parameter is essential to start chimeric read detection and reporting

Later, we identify the chimeras that align to two chromosomes, each segment longer than the threshold set with `--chimSegmentMin`. STAR-Fusion is part of the suite of tools, [Trinity Cancer Transcriptome Analysis Toolkit (CTAT)](https://github.com/NCIP/Trinity_CTAT/wiki), focused on identifying and characterizing fusion transcripts in cancer. The Trinity CTAT ecosystem of tools requires a CTAT genome lib, which is effectively a resource package containing a target genome, reference annotations, and various meta data files that support fusion-finding. On BIANCA a pre-compiled CTAT genome lib is available that I am going to use. Predicted fusions are then 'in silico validated' using FusionInspector, which performs a more refined exploration of the candidate fusion transcripts, runs Trinity to de novo assemble fusion transcripts from the RNA-Seq reads, and provides the evidence in a suitable format to facilitate visualization. Predicted fusions are annotated according to prior knowledge of fusion transcripts relevant to cancer biology (or previously observed in normal samples and less likely to be relevant to cancer), and assessed for the impact of the predicted fusion event on coding regions, indicating whether the fusion is in-frame or frame-shifted along with combinations of domains that are expected to exist in the chimeric protein.

```bash
echo "Detecting gene fusions on sample: ${sample}"
JUNCTIONS="${INPUT_DIR}${sample}-Chimeric.out.junction"
        
OUTPUT_DIR="${FUSION_DIR}${sample}/"
mkdir -p "${OUTPUT_DIR}"
STAR-Fusion \ 
  --genome_lib_dir "${CTAT}" \ 
  --chimeric_junctions "${junction_file}" \ 
  --FusionInspector validate \ 
  --denovo_reconstruct \ 
  --examine_coding_effect \ 
  --output-dir "${OUTPUT_DIR}" 
```
 Running parameters:
- `--genome_lib_dir [PATH]` – path to the CTAT genome lib,
- `-J/--chimeric_junctions [PATH]` – accession path to the 'Chimeric.out.junction' file created by STAR
- `--FusionInspector 'validate'` – requires the candidates (in a format: geneA--geneB) from the first column of the 'Chimeric.out.junction' file, then aligns them to the genome to identify those reads that align concordantly as fusion evidence to the fusion contigs,
- `--denovo_reconstruct` – *de novo* reconstruction of fusion transcripts using Trinity, creates a `.fasta` and a `.bed` file,
- `--examine_coding_effect` –  explores the effect the fusion events have on coding regions of the fused genes.

## 2.4. Create report

> [!NOTE]
> It is important to check the quality of the mapping process. Either with a reference or de novo assembly, the complete reconstruction of transcriptomes using short reads is challenging, they sometimes align equally well to multiple locations (multi-mapped reads or multi-reads). Paired-end reads reduce the problem of multi-mapping, because a pair of reads must map within a certain distance of each other and in a certain order. The percentage of mapped reads is a global indicator of the overall sequencing accuracy and of the presence of contaminating DNA.

As a last step we produce statistics to gather the necessary information about the success of the alignment with **Picard tools (v.3.1.1)**, using the following commands:
- `CollectAlignmentSummaryMetrics`- produces a summary of alignment metrics from the sorted BAM files, detailing the quality of the read alignments as well as the proportion of the reads that passed machine signal-to-noise threshold quality filters (*specific to Illumina data*),
- `CollectInsertSizeMetrics` - collect metrics about the insert size distribution of a paired-end library, useful for validating library construction including the insert size distribution and read orientation of paired-end libraries (*the expected proportions of these metrics vary depending on the type of library preparation used, resulting from technical differences between pair-end libraries and mate-pair libraries*),
- `CollectGcBiasMetrics` - collect metrics regarding GC bias: the relative proportions of guanine (G) and cytosine (C) nucleotides in a sample. Regions of high and low G + C content have been shown to interfere with mapping/aligning, ultimately leading to fragmented genome assemblies and poor coverage in a phenomenon known as 'GC bias'.

```bash
# Extract alignment metrics with Picard
${PICARD_EXE} CollectAlignmentSummaryMetrics I="${sample}" O="${picard_output}/${name}/alignment_metrics.txt"

# Extract insert size metrics with Picard
${PICARD_EXE} CollectInsertSizeMetrics I="${sample}" O="${picard_output}/${name}/insert_size_metrics.txt" \
H="${picard_output}/${name}/insert_size_histogram.pdf"

# Extract GC bias metrics with Picard
${PICARD_EXE} CollectGcBiasMetrics I=${sample} O="${picard_output}/${name}/gc_bias_metrics.txt" \
CHART="${picard_output}/${name}/gc_bias_chart.pdf" S="${picard_output}/${name}/gc_summary.txt"
```

Last but not least, we parse summary statistics from results and log files generated by all other bioinformatics tools with **MultiQC (v.1.10.1)**. MultiQC recursively searches through any provided file paths and finds files that it recognises. It parses relevant information from these and generates a single stand-alone HTML report file. Furthermore, a multiqc_data folder will be generated as well with the reformatted compact data tables in there, one from each module in TSV format, as well as a verbose multiqc.log file and a multiqc_data.json. MultiQC's highly collaborative modules can work with the log files produced during the analysis:
- **Cutadapt** – this module summarizes found and removed adapter sequences, primers, poly-A tails and other types of unwanted sequences, 
- **FastQC** – the quality control module generates similar plots, which are included in the default HTML report as well. An additional fastqc_data.txt is generated too that can be helpful for downstream analysis as it is relatively easy to parse,
- **STAR** - module parses summary statistics from the *Log.final.out* log files,
- **featureCounts** – module parses results generated by featureCounts, visualizes the reads mapped to genes, exons, promoter, gene bodies, genomic bins, chromosomal locations or other features. The filenames must end in *'.summary'* to be discovered,
- **Hisat2** – module parses summary statistics if option `--new-summary` has been specified,
- **Salmon** – module parses *meta_info.json, lib_format_counts.json* and *flenDist.txt* files, if found.

Whilst MultiQC is typically used as a final reporting step in an analysis, it can also be used as an intermediate in your analysis, as these files essentially standardize the outputs from a lot of different tools and make them useful for downstream analysis. We will do that as well, going one step further with the TidyMultiqc package. The TidyMultiqc package provides the means to convert the multiqc_data.json file into a tidy data frame for downstream analysis in R. 