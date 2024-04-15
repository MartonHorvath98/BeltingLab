# 1. Integrative mRNA and Proteomic Analysis of Glioblastoma cell-surface antigens

$${\color{gray}Project\ works\ done\ for\ Mattias\ Belting's\ laboratory\ at\ Lund\ University}$$

- [1. Integrative mRNA and Proteomic Analysis of Glioblastoma cell-surface antigens](#1-integrative-mrna-and-proteomic-analysis-of-glioblastoma-cell-surface-antigens)
  - [1.1. Introduction](#11-introduction)
  - [1.2. Project aims](#12-project-aims)
  - [1.3. Sample Metadata Summary](#13-sample-metadata-summary)
    - [1.3.1. Metadata Fields](#131-metadata-fields)
  - [1.4. Setting Up the Environment](#14-setting-up-the-environment)
    - [1.4.1. Directory tree](#141-directory-tree)
    - [1.4.2. Environment prerequisits](#142-environment-prerequisits)
  - [1.5. Raw data files](#15-raw-data-files)
- [2. Bioinformatical pipeline](#2-bioinformatical-pipeline)
  - [2.1. Quality Control](#21-quality-control)
  - [2.2. Prepare reference genome](#22-prepare-reference-genome)
  - [2.3. Read mapping](#23-read-mapping)
    - [2.3.1 HISAT2](#231-hisat2)
    - [2.3.2. Salmon](#232-salmon)
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
### 1.4.2. Environment prerequisits

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

## 1.5. Raw data files

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

# 2. Bioinformatical pipeline 
$${\color{gray}-\ Compare\ gene\ expression\ RNA-seq\ from\ cell\ lines\ (2D)\ and\ organoids\ (3D)\ grown\ under\ normoxia\ and\ hypoxia\ -}$$

## 2.1. Quality Control
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

## 2.2. Prepare reference genome

Once the pre-processing and quality control steps are completed the resulting, high-quality data is ready for the read mapping or alignment step. Depending on the availability of a reference genome sequence, it can happen in one of two ways: 
1. When studying an organism with a reference genome, it is possible to infer which transcripts are expressed by mapping the reads to the reference genome (Genome mapping) or transcriptome (Transcriptome mapping). Mapping reads to the genome requires no knowledge of the set of transcribed regions or the way in which exons are spliced together. This approach allows the discovery of new, unannotated transcripts.
2. When working on an organism without a reference genome, reads need to be assembled first into longer contigs (de novo assembly). These contigs can then be considered as the expressed transcriptome to which reads are re-mapped for quantification.

>[!NOTE]
>*In this case, we follow the aprroach remains genome mapping, followed by building a more accurate, splicing-sensitive transcriptome in a genome-guided manner, and then to increase the number of mapped reads by aligning them to the assembled transcriptome. However, these steps will only happen as part of the later workflow to identify gene fusions, and mutated genes.*

There are many bioinformatics tools available to perform the alignment of short reads. Here, the one used is **Hisat2 (v.2.2.1)**, one of the most popular mapper tools, which stands for “hierarchical indexing for spliced alignment of transcripts 2”. Hisat2 provides fast and sensitive alignment, in a splice-sensitive manner, which makes it ideal for aligning RNA-samples. Besides Hisat2, read will also be aligned with quasi-mapper, **Salmon (v.1.10.1)**. The technique employed by Salmon rapidly approximates the positions where reads might originate from the transcriptome without producing an explicit alignment of reads to the genome allowing a rapid quantification of gene expression levels.

>[!NOTE]
>*Besides the mapper algorithms, we will use Samtools for processing and analysing the data. It includes tools for file format conversion and manipulation, sorting, querying, and statistics amongst other methods.*

Indexes for the hisat2 and salmon aligners are created using the `makeGrch38.sh` script. The hisat2-build command will create eight files with the extensions *.1.ht2, .2.ht2, .3.ht2, .4.ht2, .5.ht2, .6.ht2, .7.ht2* and *.8.ht2*. These files together are the index for hisat2. Because the plans include heavy focus on alternative splicing and gene fusion events, the genome index is created with transcript annotations using the `--ss` and `--exon` flags. These take into account the prepared exons and splice sites, in HISAT2's own format (three or four tab-separated columns).

>[!IMPORTANT]
>***It is crucial to download the same reference in .fasta and .gff format from the same origin, to avoid conflicts in later steps of the analysis pipeline. Here, the GRCh38 release version v.108 (2022 Oct) of H. sapiens (human) from Ensembl is used as reference.***

```bash
# The Ensembl release version and base URLs for downloads
ENSEMBL_RELEASE=108
ENSEMBL_GENOME=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna
ENSEMBL_GFF3_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gff3/homo_sapiens/
ENSEMBL_TRANSCRIPTOME=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/cds
# Index reference genome with samtools 
samtools faidx <path/…/ref.fasta> 
# Convert GFF file to GTF 
gffread <path/…/ref.gff> -T -o <path/…/ref.gtf>
# Extract splice sites and exons
hisat2_extract_splice_sites.py -v <path/…/ref.gtf> > <path/…/splicesite.tsv>
hisat2_extract_exons.py -v <path/…/ref.gtf> > <path/…/exons.tsv>

# Build HISAT2 index with transcripts annotations
hisat-build -p ${threads} --seed 100 ${ENSEMBL_GENOME} --ss <path/…/splicesite.tsv> --exon  <path/…/exons.tsv>  <path/…/HISAT/base_name>
```

To get accurate quantification estimates with Salmon in mapping-based mode a salmon index for the transcriptome has to be built as well. We build the recommended, decoy-aware transcriptome file, using the entire genome of the organism (Grch38) as the decoy sequence. This is achieved by concatenating the genome to the end of the transcriptome before indexing. Salmon indexing also requires populating the decoys.txt file with the chromosome names, which is extractable by using the `grep` command. This scheme provides a more comprehensive set of decoys, but, obviously, requires considerably more memory to build the index. This will build the mapping-based index, using an auxiliary k-mer hash (*default: k-mers of length 31*), that will also act as the minimum acceptable length for a valid match. Running an exploratory mapping step, we found low alignment rate (= appx. 40%), thus a smaller value of k is set to improve sensitivity even more using selective alignment in the later steps (enabled via the `–validateMappings` flag). We can run salmon indexing step, with the same ensemble reference files, as follows:

>[!IMPORTANT] 
> The genome targets (decoys) should come after the transcriptome targets in the reference!

```bash
# Create decoys.txt file
grep "^>" ${ENSEMBL_GENOME} | cut -d " " -f 1 | sed -e 's/>//g' > <path/…/decoys.txt>
# Concatenate transcriptome and genome
cat ${ENSEMBL_TRANSCRIPTOME} ${ENSEMBL_GENOME} > <path/…/gentrome.fa>

# Build Salmon index
salmon index -t gentrome.fa -d decoys.txt -i <path/…/SALMON/> -k 13 -p ${threads}
```

## 2.3. Read mapping

> [!IMPORTANT]
> **The most important parameter to get right is the strandedness of the library. Salmon can predict the library preparation method when used with the automatic library preparation tag `--libType "A"`. Executing an exploratory run with Salmon told that the library in unstranded (*"IU"*), however there is a very high mapping bias towards the reverse as in the case of more than 99% of the individual alignments the first read maps to the reverse strand (*"ISR"*). Consequentially, during the hisat2 alignment we use the option `--rna-strandedness "RF"`, that corresponds to the second-strand synthesis method (common for many Illumina kits), where R1 maps to the reverse strand of the DNA. During the Salmon qunatification we sat --libType 'ISR'.**

### 2.3.1 HISAT2

We align the RNA-seq reads to a genome in order to identify exon-exon splice junctions with the help of the `alignReads.sh` script. The next downstream step would most likely be a genome-guided transcriptome assembly, so we use the --novel-splicesite-outfile mode to reports a list of splice sites in the file: chromosome name tab genomic position of the flanking base on the left side of an intron tab genomic position of the flanking base on the right tab strand (+, -, and .) ‘.’ indicates an unknown strand for non-canonical splice sites. During the alignment step the applied parameters mean:
- `--phred33` – defines the input quality (*default: Phred+33*),
- `--dta` – (=downstream-transcriptome-assembly) report alignments tailored for transcript assemblers including StringTie,
- `--p [INT]` – sets INT number of cores to be used during the calculations,
- `--non-deterministic` – Hisat2 re-initializes its pseudo-random generator for each read using the current time, meaning that Hisat2 will not necessarily report the same alignment for two identical reads. This might be more appropriate in situations where the input consists of many identical reads,
- `--rna-strandedness 'RF'` – is used to define the strand-specificity of the RNA-seq protocol,
- `--x` – the path to the index for the reference genome including the basename of the files excluding the *.1.ht2 / etc.* extension,
- `--1` and `--2` – are comma-separated list of files containing mate 1s and mate 2s,
- `--summary-file` – sets the path to the summary file,
- `--new-summary` – makes the summary file readable by MultiQC for downstream analysis.

Following the alingment step, **Samtools (v.1.19.2)**, channeled through the piping operator '|', instantly converts the .sam files produced by the aligner to their binary counterpart .bam, which limits space requirements as the output file size shrinks considerably via `samtools view -h -bS > <path/…/file.bam>`. Then using the `samtools sort` and `index` command the bam files are sorted along the reference, and BAI-format index is generated. Finally,  **featureCounts (v.2.0.6)** takes as input the BAM files and the genome annotation file in GFF format containing the chromosomal coordinates of features. It outputs numbers of reads assigned to features or meta-features. Each entry in the GFF annotation file is considered a feature (e.g. an exon), while a meta-feature is the aggregation of a set of features (e.g. a gene). When summarising reads at meta-feature level, read counts obtained for features included in the same meta-feature will be added up to yield the read count for the corresponding meta-feature. It also outputs stat info for the overall summarization results, including number of successfully assigned reads and number of reads that failed to be assigned due to various reasons (these reasons are included in the stat info). While counting how many reads align to each gene in a genome annotation the set parameters are:
- `-p` – specifies that input data contain paired-end reads and performs fragment counting (ie. counting read pairs), the `--countReadPairs` parameter should also be specified in addition to this parameter
- `-s [INT]` – determines the strand-specificity of the read counting (*here: 2 (reversely stranded)*),
- `-f` – performs read counting at feature level (eg. counting reads for exons rather than genes),
- `-O` –  assigns reads to all their overlapping features (or meta-features if -f is not specified),
- `-M` – multi-mapping reads will also be counted,
- `-T [INT]` – specifies the number of threads to be used,
- `-a` – is a mandatory argument, the path to the annotation file,
- `-o` –  is the path to and name of the output file that includes the read counts and a separate file with summary statistics.

```bash
 # Align reads with Hisat2
hisat2 --phred33 --dta --non-deterministic -p ${cores} --novel-splicesite-outfile "${output_dir}/${sample}/novel_splicesite.txt"\
--summary-file "${output_dir}/${sample}/stats.txt" --new-summary -x ${reference}/hisat/${base_name} -1 ${fw_read} -2 ${rv_read} |\
samtools view -h -bS > "${output_dir}/${sample}.bam"

# Sort the bam file
samtools sort -@ ${cores} -o "${output_dir}/${sample}.sorted.bam" "${output_dir}/${sample}.bam"

# Index the sorted bam file
samtools index "${output_dir}/${sample}.sorted.bam"

# Remove the unsorted bam file
rm "${output_dir}/${sample}.bam"

# Calculate readcounts with featureCounts
featureCounts -p --countReadPairs -s 0 -f -M -O -T ${cores} -a ${reference}/Homo_sapiens.GRCh38.gff\
-o "${output_dir}/${sample}.counts.txt" "${output_dir}/${sample}.sorted.bam"
```

### 2.3.2. Salmon

We quantify the reads directly in mapping-based mode against the prepared index using the `salmon quant` command. We introduced selective alignment with the `--validateMappings` flag to adopt a considerably more sensitive scheme for finding the potential mapping loci of a read, and score potential mapping loci using the dynamic programming chaining algorithm, SIMD, introduced in minimap2. We emply the pre-pepared decoy-aware transcriptome index to mitigate potential spurious mappings to unannotated genomic locuses with high sequence-similarity to the annotated transcriptome. For the mapping-based quantification the set parameters are:
- `-l 'ISR'` – sets the library type to inward orientation ('I'), stranded ('S'), reverse-first ('R') protocol,
- `-i [PATH]` –  the path to the reference genome index,
- `-1` and `-2` –  mate 1s and mate 2s reads,
- `-p [INT]` – specifies the number of threads to be used,
- `--seqBias` – performs sequence-specific bias correction,
- `--validateMappings` – calls selective alignment,
- `-o` – is the path to the output folder where the quantification file (called quant.sf), the command information file (called cmd_info.json) and a file (called lib_format_counts.json) with the number of fragments that had at least one mapping compatible with the designated library format are saved.

```bash
# Align reads with Salmon
salmon quant -l 'ISR' -i ${reference} -1 ${fw_read} -2 ${rv_read} -p ${threads} --seqBias --validateMappings -o ${output_dir}/${sample}
```

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
- **Hisat2** – module parses summary statistics if option `--new-summary` has been specified,
- **featureCounts** – module parses results generated by featureCounts, visualizes the reads mapped to genes, exons, promoter, gene bodies, genomic bins, chromosomal locations or other features. The filenames must end in *'.summary'* to be discovered.

Whilst MultiQC is typically used as a final reporting step in an analysis, it can also be used as an intermediate in your analysis, as these files essentially standardize the outputs from a lot of different tools and make them useful for downstream analysis. We will do that as well, going one step further with the TidyMultiqc package. The TidyMultiqc package provides the means to convert the multiqc_data.json file into a tidy data frame for downstream analysis in R. 