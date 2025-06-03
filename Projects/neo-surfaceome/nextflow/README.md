# Bioinformatic Pipeline for RNA-seq data processing

- [Bioinformatic Pipeline for RNA-seq data aprocessing](#bioinformatic-pipeline-for-rna-seq-data-processing)
  - [1. Why nextflow?](#1-why-nextflow)
  - [2. Nextflow Configuration](#2-nextflow-configuration)
    - [2.1. Profiles](#21-profiles)
    - [2.2. Processes](#22-processes)  
    - [2.3. Parameters](#23-parameters)
  - [3. Subworkflows](#3-subworkflows)
    - [3.1. Core Workflow](#31-core-workflow)
    - [3.2. Alignment Workflow](#32-alingnment-workflow)
    - [3.3. Fusion Workflow](#33-fusion-workflow)
    - [3.4. Variant Workflow](#34-variant-workflow)
  - [4. Complete workflow](#4-complete-workflow)
    - [4.1. Running the Pipeline](#41-running-the-pipeline)

# 1. Why nextflow?

For this project we used **Nextflow** *(v24.10.6)* [(Di Tomasso P et al., 2017)](https://doi.org/10.1038/nbt.3820) to ensure reproducibility and scalability. Nextflow is optimized for parallelization and distributed computing, that made it an optimal tool for us to tackle our *'lower than required'* throughput problem. For installation, environment setup and developement guidelines the official documentation of [Nextflow](https://www.nextflow.io/docs/latest/index.html) provided an invaluable resource.

Shortly, Nextflow is a workflow language, that runs on the Java virtual machine (JVM), it has a syntax very similar to Groovy, allowing for making full use of the Java and Groovy standard libraries. Nextflow follows the Unix philosophy, in which many simple command line tools can be chained together into increasingly complex tasks. Nextflow pipeline script, can include different processes written in any scripting language (Bash, Perl, Ruby, Python, etc.), which was very beneficial for us, so that our previous scripts could have been integrated with minimal updates. In practice, these processes are then executed independently, and communicate only via asynchronous FIFO channels (=queues), from which we can define one or more as input and output.

To set up an optimal development environment ofr creating, testing, and optimizing the Nextflow pipeline we used [VS Code](https://code.visualstudio.com/download), with the VS CODE [Nextflow extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow), to add language support and syntax highlights enhancing, diagnostics and debugging.

> [!IMPORTANT]
> ***Nextflow uses UTF-8 as the default character encoding for source files, we have to make sure to use UTF-8 encoding when editing Nextflow scripts with VS-code and keep it in mind as a possible debug option!***

# 2. Nextflow Configuration

The `nextflow.config` file is a central configuration file used by Nextflow to customize and control the behavior of workflows. It allows us to define parameters, execution settings, environment profiles, and resource requirements without hardcoding them into the pipeline scripts. When a pipeline script is launched, Nextflow looks for configuration files in multiple locations (by default), and apply them from lowest to highest priority (to avoid conflicting settings). In this project, however, we generally used the `-C <config-file>` command line option to specify a fixed configuration file `$NXF_HOME/nextflow.config` and ignore all other default file locations.

```bash
nextflow run complete-workflow.nf -C $NXF_HOME/nextflow.config
```

## 2.1. Profiles

Profiles in Nextflow are named sets of configuration settings that let us easily switch between different computing environments, resource configurations, or modes of running the pipeline.

1. **'local' profile — for development and testing:**

A local profile was defined for testing the pipeline on my personal computer: it uses the local executor - 'bash' on my virtual machine (WSL2) and minimal resources:

```groovy
local {
    docker.enabled = true

    process {
        executor = 'local'
    }
}
```

2. **'uppmax' profile - for production data processing:**

The uppmax profile is configured to run on the [BIANCA](https://www.uu.se/en/centre/uppmax/resources/clusters/bianca) HPC cluster, reserved for NAISS-SENS projects, i.e. sensitive data, on UPPMAX (Uppsala Multidisciplinary Center for Advanced Computational Science). UPPMAX is a shared resource, and it uses the [SLURM](https://docs.uppmax.uu.se/cluster_guides/slurm/) sceduling system to ensure fair allocation. SLURM can be accessed either by submitting and running a batch job script via the `sbatch` command or by starting an interactive session with `interactive -A [project-code]`. The sbatch script must include Slurm directives and the commands necessary to execute the job:
- **mandatory settings**: 
  - *account (-A)*: you must be a member to charge a project 
- **important settings**:  
  - *partition (-p)*: should the project run of complete nodes or only parts of them, on cores (default!) (1 node = 16 core(s)). The core hours scale with the number of course, which is important to keep in mind not to overshoot the 5000 core-hours allocated to a project each month.
  - *number of cores/nodes (-n)*
  - *time (-t)*: estimated runtime in core hours (hh:mm:ss). 

> [!IMPORTANT]
> ***Always overestimate with ~50%, because jobs get killed when they reach the time limit, but users only get charged for the actual runtime.*** 

Furthermore, on the shared Linux computers of the UPPMAX clusters users cannot modify, upgrade or uninstall software themselves and instead an [Environment Module System](https://lmod.readthedocs.io/en/latest/) is used. Luckily Nextflow provied module system support, but it has to be set up in th econfig file as well:
```groovy
uppmax {
        params.max_memory = 500.GB          //defaults for BIANCA
        params.max_cpus = 16                //defaults for BIANCA
        params.max_time = 240.h             //defaults for BIANCA
	    params.project = 'sens1234123'      //Replace with actual project ID
	    executor.account = 'sens1234123'    //Replace with actual project ID
	    executor.queueSize = 60

        process {
            executor = 'slurm'
            scratch = '$SNIC_TMP'
        }
}
```
Nextflow provides an abstraction between the pipeline's functional logic and the underlying execution system, thus the same process can get the next two example headers based on the profile selected on the command line from the configuration file:
```bash
nextflow run complete-workflow -C $NXF_HOME/nextflow.config -profile 'local'
```
The example header:
```bash
#!/bin/bash -l
```
BUT, if:
```bash
nextflow run complete-workflow -C $NXF_HOME/nextflow.config -profile 'uppmax'
```
The example header:
```bash
#!/bin/bash -l
#SBATCH –A sens1234123 
#SBATCH –J tag 
#SBATCH –p core 
#SBATCH –n 4 
#SBATCH –t 06:00:00 
```

## 2.2 Processes

In Nextflow, a process is the fundamental unit of computation. Each process defines a self-contained task (*e.g. running a script*), and is executed independently according to the data dependencies defined by the workflow. Processes are portable, parallelizable, and reproducible, making them ideal for building scalable workflows across local machines, HPC clusters, and cloud platforms.

A process consists of a name and a body. The process body define additional sections for directives, inputs, outputs, **script** (compulsory!), etc. 

> [!NOTE]
> A process must define a script block, all other sections are optional. Directives do not have an explicit section label, but must be defined first (see below). For example `modules/fastqc/main.nf`:

```groovy
process FASTQC {
    tag "FastQC on $sample" //#SBATCH –J tag
    label "low_effort" // #SBATCH –p core, –n 4, –t 06:00:00 

    publishDir "${params.report_outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)
    val subfolder // to avoid overwriting pre-trim with post-trim QC

    output:
    tuple val(sample), path("$subfolder/${sample}/*_fastqc.zip") , emit: fastqc_report
    path "versions.yml"                                          , emit: version

    script:
    """
    mkdir -p $subfolder/${sample}
    fastqc \\
        -o $subfolder/${sample} \\
        ${params.fastqc_args} \\
        --threads ${task.cpus} \\
        ${read1} ${read2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}
```

In the two profiles 'local' and 'uppmax' based on the directives, two different set of process selectors were set up to define executor specific properties for every process. The `withLabel` selectors allow the configuration of all processes annotated with a label directive, while the `withName` selector allows the configuration of a specific process in the pipeline by its name. In our pipeline, labels defined the time, memory and/or wall-clock time allocation for each process, while the name selectors ensured loading the suitable environments via [Docker](https://www.docker.com/products/docker-desktop/) (v28.0.4) containers ([BioContainers](https://biocontainers.pro/)) or Environment Modules:
```groovy
profiles {
    // Local settings - executor: WSL2, bash
    local {
        process {
            withLabel: low_effort {
                cpus = 4
                memory = 4.GB
            }
            withLabel: generic_task {
                cpus = 4
                memory = 12.GB
            }
            withLabel: thicc {
                cpus = 4
                memory = 22.GB
            }
            withName:FASTQC {
                container = 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
            }
            // Many more tools...
        }
    }
    // UPPMAX settings - executor: slurm
    uppmax {
        process {
            withLabel: low_effort {
                cpus = 4
                time = { 4.h * task.attempt } //Double time in case of server timeout
                errorStrategy = { task.exitStatus == 140 ? 'retry' : 'terminate' } 
                maxRetries = 2
            }
            withLabel: generic_task {
                cpus = 8
                // Same strategy for time, etc.
            }
            withLabel: thicc {
                cpus = 16
                // Same strategy for time, etc.
            }
            withName:FASTQC {
                module = 'bioinfo-tools:FastQC/0.11.9'
            }
            // Many more tools...
        }
    }
}
```
> [!IMPORTANT]
> To access bioinformatics tools on UPPMAXX, the bioinfo-tools module must be loaded first!

## 2.3 Parameters

We define pipeline paramteres in the config file using the `params` scope. Users can overwrite any parameter while running the workflow on the command line.

> [!NOTE]
> In this project, the parameters are predifined to use the CTAT-genome-lib reference and the files included in it by default. Any other reference must be specified on the command line. Tool runtime paramteres are also specified in advance.

> [!IMPORTANT]
>  Parameters can be specified on the command line using double-dash (`--`). Single dash (`-`) parameters are used for nextflow interanlly. Otherwise, multiple parameters can also be redifined using a params file using the `-params-file` option.

# 3. Subworkflows 

[:rewind:](../README.md#bioinformatical-pipeline) *Return to the main README file, if you used the link to the Nextflow configurations, before proceeding...*  

## 3.1. Core Analysis

This subworkflow integrates (i) pre-trim QC using FastQC (v0.12.1), (ii) trimming using the TrimGalore suite (v0.6.10), (iii) post-trim QC with FastQC and SeqKit (v2.10), and (iv) isoform-level quantification using Salmon (v1.10.1) in mapping-based mode.

### Quality check (QC) - `modules/fastqc/main.nf`

```groovy
process FASTQC {
    script:
    """
    mkdir -p $subfolder/${sample}
    fastqc \\
        -o $subfolder/${sample} \\
        ${params.fastqc_args} \\
        --threads ${task.cpus} \\
        ${read1} ${read2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}
```
Transcribed to:
```bash
# e.g. pre-trim on VI-3429-593-2DH-1
mkdir -p pre-trim/VI-3429-593-2DH-1

trim_galore \
    -O pre-trim/VI-3429-593-2DH-1 \
    -q 20 --length 36 --paired --illumina \
    VI-3429-593-2DH-1_R1_001.fastq.gz VI-3429-593-2DH-1_R2_001.fastq.gz
```
*e.g.: FastQC over-represented sequences fails VI-3429-593-2DH-1_R1_001.fastq.gz:*
| #Sequence                                          | Count | Percentage | Possible Source                          |
|----------------------------------------------------|-------|------------|------------------------------------------|
| GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGGAGCATCTCGTAT | 1058610   | 1.07       | TruSeq Adapter, Index 18 (97% over 38bp) |
| GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGGAGCATCGCGTAT | 389043    | 0.39       | TruSeq Adapter, Index 18 (97% over 38bp) |
| GTCGGCATGTATTAGCTCTAGAATTACCACAGTTATCCAAGTAGGAGAGG | 165859    | 0.17       | No Hit                                  |

Following a prelimnary quality assessment with FastQC, we can say that the overall quality of the bases is high, over the required treshold, however there is a high precentage of adapter contamination. Every adapter match seems to fall under the Illumina adapters' list, so the flag `--illumina` will be included in the trimming! 

### Trimming - modules/trim_gelore/main.nf


Quality assessment is performed using the TrimGalore suite (v0.6.1), which is a wrapper script around Cutadapt, a semi-global aligner algorithm (also called free-shift). It means that the sequences are allowed to freely shift relative to each other and differences are only penalised in the overlapping region between them. The algorithm works using unit costs (alignment score) to find the optimal overlap alignments, where positive value is assigned to matching bases and penalties are given for mismatches, inserts or deletions.

>[!IMPORTANT]
> ***It is important to check that sequence quality is similar for all samples and discard outliers. As a general rule, read quality decreases towards the 3’ end of reads, and if it becomes too low, bases should be removed to improve mappability. The quality and/or adapter trimming may result in very short sequences (sometimes as short as 0 bp), and since alignment programs may require sequences with a certain minimum length to avoid crashes to short fragments (in the case above, below 36 bases: --length 36) should not be considered either.***


### Isoform-level quantification

For read quantification we used the Salmon (v1.10.1), a fast and lightweight method for quantifying transcript abundance from RNA–seq reads. 
Prepare reference genome



## 2.3. Read mapping

>[!NOTE]
>*Besides the STAR mapper algorithm, Samtools is used for processing and analysing the data. It includes tools for file format conversion and manipulation, sorting, querying, and statistics amongst other methods.*

We align the RNA-seq reads to a genome with STAR in order to identify abnormal splicing events and chimeras (fusion genes) with the help of the `sbatch_STAR.sh` script. Following the alignment we will use **STAR-fusion (v1.10.1)** to identify candidate fusion transcripts supported by the Illumina reads.

```bash
# Read the trimmed forward and reverse reads
echo "Aligning reads with STAR on sample: ${sample}"
    
fw_read=${INPUT_DIR}${sample}_R1_001_val_1.fq.gz
rv_read=${INPUT_DIR}${sample}_R2_001_val_2.fq.gz

# Align reads with STAR
STAR --runThreadN ${CORES}
  --twopassMode Basic
  --genomeDir ${REF} 
    --readFilesIn ${fw_read} ${rv_read}
    --readFilesCommand "gunzip -c"
    --outSAMtype BAM SortedByCoordinate
    --outBAMcompression 6
    --outSAMstrandField intronMotif
    --outSAMunmapped Within
    --outFileNamePrefix "${OUTPUT_DIR}${sample}-"
    --outReadsUnmapped None
    --outFilterScoreMinOverLread 0
    --outFilterMatchNminOverLread 0
    --outFilterMatchNmin 0
    --outFilterMismatchNmax 2
    --alignSJstitchMismatchNmax 5 -1 5 5
    --alignIntronMin 10
    --alignIntronMax 100000
    --alignMatesGapMax 100000
    --chimOutType Junctions # **essential** to create the Chimeric.junction.out file for STAR-Fusion
    --chimSegmentMin 12
    --chimJunctionOverhangMin 8 # **essential to invoke chimeric read detection & reporting**
    --chimOutJunctionFormat 1 # **essential** includes required metadata in Chimeric.out.junction file.
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

## 2.4. Fusion detection

In the next step of analysis we identify the chimeras that align to two chromosomes, each segment longer than the threshold set with `--chimSegmentMin`. STAR-Fusion is part of the suite of tools, [Trinity Cancer Transcriptome Analysis Toolkit (CTAT)](https://github.com/NCIP/Trinity_CTAT/wiki), focused on identifying and characterizing fusion transcripts in cancer. The Trinity CTAT ecosystem of tools requires a CTAT genome lib, which is effectively a resource package containing a target genome, reference annotations, and various meta data files that support fusion-finding. On BIANCA a pre-compiled CTAT genome lib is available that I am going to use. Predicted fusions are then 'in silico validated' using FusionInspector, which performs a more refined exploration of the candidate fusion transcripts, runs Trinity to de novo assemble fusion transcripts from the RNA-Seq reads, and provides the evidence in a suitable format to facilitate visualization. Predicted fusions are annotated according to prior knowledge of fusion transcripts relevant to cancer biology (or previously observed in normal samples and less likely to be relevant to cancer), and assessed for the impact of the predicted fusion event on coding regions, indicating whether the fusion is in-frame or frame-shifted along with combinations of domains that are expected to exist in the chimeric protein.

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
- **featureCounts** – module parses results generated by featureCounts, visualizes the reads mapped to genes, exons, promoter, gene bodies, genomic bins, chromosomal locations or other features. The filenames must end in *'.summary'* to be discovered.

Whilst MultiQC is typically used as a final reporting step in an analysis, it can also be used as an intermediate in your analysis, as these files essentially standardize the outputs from a lot of different tools and make them useful for downstream analysis. We will do that as well, going one step further with the TidyMultiqc package. The TidyMultiqc package provides the means to convert the multiqc_data.json file into a tidy data frame for downstream analysis in R. 

[!IMPORTANT]
>***When logging in, all users join the login node first. It is IMPORTANT, although not compulsory on BIANCA, to move onto a computing node (using either an interactive session or an sbatch script) to perform calculations. However, memory intensive processes will be down prioritized to not block access for others (and might never be finished).***