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
nextflow run main.nf -C $NXF_HOME/nextflow.config -profile 'local'
```
The example header:
```bash
#!/bin/bash -l
```
BUT, if:
```bash
nextflow run main.nf -C $NXF_HOME/nextflow.config -profile 'uppmax'
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

:rewind: *[Return](../README.md#bioinformatical-pipeline) to the main README file, if you used the link to the Nextflow configurations, before proceeding...* 

Subworkflows are named workflows that can be called by other workflows, in our case only by the entry workflow (`complete-workflow.nf`). In named (sub)workflows a `take:` section is used to declare the input channels, and the `emit:` section is used to declare the output channels. Inputs can be specified just like arguments when calling the workflow and the output channels can be accessed using the `out` property similarly to how it works in processes.

## 3.1. Core workflow

This subworkflow integrates:
1. pre-trim QC using FastQC (v0.12.1) and SeqKit (v2.10)
2. trimming using the TrimGalore suite (v0.6.10) 
3. post-trim QC with FastQC and SeqKit
4. isoform-level quantification using Salmon (v1.10.1) in mapping-based mode

### Quality check (QC) - `modules/fastqc/main.nf`

We used [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (v0.12.1) to execute simple quality control checks on raw and trimmed sequence data. It provides a modular set of analyses: importing FastQ files, providing a quick overview to tell us about problems, and creating summary graphs and tables to quickly assess the data.

```groovy
process FASTQC {
    // Directives
    // input:
    // output:
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
# e.g. raw_data on VI-3429-593-2DH-1
mkdir -p raw_data/VI-3429-593-2DH-1

fastqc \
    -O raw_data/VI-3429-593-2DH-1 \
    -f fastq -q \
    VI-3429-593-2DH-1_R1_001.fastq.gz VI-3429-593-2DH-1_R2_001.fastq.gz
```
*e.g.: FastQC over-represented sequences fails VI-3429-593-2DH-1_R1_001.fastq.gz:*
| #Sequence                                          | Count | Percentage | Possible Source                          |
|----------------------------------------------------|-------|------------|------------------------------------------|
| GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGGAGCATCTCGTAT | 1058610   | 1.07       | TruSeq Adapter, Index 18 (97% over 38bp) |
| GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCGGAGCATCGCGTAT | 389043    | 0.39       | TruSeq Adapter, Index 18 (97% over 38bp) |
| GTCGGCATGTATTAGCTCTAGAATTACCACAGTTATCCAAGTAGGAGAGG | 165859    | 0.17       | No Hit                                  |

Following a prelimnary quality assessment with FastQC, we can say that the overall quality of the bases is high, over the required treshold, however there is a high precentage of adapter contamination. Every adapter match seems to fall under the Illumina adapters' list, so the flag `--illumina` will be included in the trimming! 

### Trimming - `modules/trim_galore/main.nf`

Quality assessment is performed using the [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) suite (v0.6.1), which is a wrapper script around Cutadapt, a semi-global aligner algorithm (also called free-shift). It means that the sequences are allowed to freely shift relative to each other and differences are only penalised in the overlapping region between them. The algorithm works using unit costs (alignment score) to find the optimal overlap alignments, where positive value is assigned to matching bases and penalties are given for mismatches, inserts or deletions.

```groovy
process TRIM_GALORE {
    // Directives
    // input:
    // output:
    script:
    """    
    trim_galore ${read1} ${read2} \\
     -j $task.cpus \\
     --gzip \\
     $params.trim_args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
```
Transcribed to:
```bash
# e.g. VI-3429-593-2DH-1
trim_galore VI-3429-593-2DH-1_R1_001.fastq.gz VI-3429-593-2DH-1_R2_001.fastq.gz\
    -J 4 \
    --gzip \
    -q 20 --length 36 --paired --illumina
```
*e.g.: Trim galore log of sample VI-3429-593-2DH-1, trimming  \*R1_001.fastq.gz & \*R2_001.fastq.gz:*
| === Summary === | |
| --- | --- |
| Total reads processed: | 99,262,678 |
| Reads with adapters: | 61,605,480 (62.1%) |
| Reads written (passing filters): | 99,262,678 (100.0%) |
| Total basepairs processed: | 14,988,664,378 bp | 
| Quality-trimmed: | 19,015,148 bp (0.1%) | 
| Total written (filtered): | 13,252,226,643 bp (88.4%) |

>[!IMPORTANT]
> ***It is important to check that sequence quality is similar for all samples and discard outliers. As a general rule, read quality decreases towards the 3’ end of reads, and if it becomes too low, bases should be removed to improve mappability. The quality and/or adapter trimming may result in very short sequences (sometimes as short as 0 bp), and since alignment programs may require sequences with a certain minimum length to avoid crashes to short fragments (in the case above, below 36 bases: --length 36) should not be considered either.***


### Isoform-level quantification via Salmon - `module/salmon/quant/main.nf`

For read quantification we used the [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) (v1.10.1), a fast and lightweight method for quantifying transcript abundance from RNA–seq reads. 

```groovy
process SALMON_QUANT {    
    // Directives
    // input:
    // output:
    script:
    """
    salmon quant --threads $task.cpus \\
        $params.quant_args \\
        -i $salmon_index \\
        -1 ${trim1} -2 ${trim2} \\
        -o $sample_id

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
```
Transcribed to:
```bash
# e.g. VI-3429-593-2DH-1
salmon quant --threads 16 \\
        -l ISR \\
        -i ref_genome.salmon.idx \\
        -1 VI-3429-593-2DH-1_R1_001_val_1.fq.gz -2 VI-3429-593-2DH-1_R2_002_val_2.fq.gz \\
        -o VI-3429-593-2DH-1
```
> [!IMPORTANT]
> It is important to learn the strandedness of the experiment. Here - `ISR` - we knew that the used protocol, Illumina, preserve the original RNA sequence, that is the first strand. If we are not aware of the method and the strandedness used during the sequencing, before starting the data analysis we can learn about the directionality of the reads using the `infer_experiment.py` algorithm from the [RSeQC](https://rseqc.sourceforge.net/) program.

> [!NOTE]
> Actually, the Salmon developer team recommends allocating only 8 — 12 threads, to reach the maximum speed, threads allocated above this limit will likely spend most of their time idle / sleeping.

#### Preparing a genome index - `module/salmon/index/main.nf`

To run Salmon in mapping-based mode, first we have to build a salmon index using a decoy-aware transcriptome file, that we can create by creating a decoy file:
```bash
grep "^>" ${genome} | cut -d " " -f 1 | sed -e 's/>//g' > decoys.txt    
``` 
And concatenating the transcriptome and genome files, to run the index building on this generated `gentrome.fa`. The internal function `get_index_channel()` checks if the index has previously been created, otherwise calls appropriate process:
```groovy 
def get_index_channel() {
    def index_path = file(params.salmon_index)
    return index_path.exists()
        ? Channel.value(index_path)
        : INDEX(params.transcriptome_file, params.genome_file) }
```

## 3.2. Alignment workflow

This subworkflow integrates:
1. splice-aware genome-guided alignment using STAR (v2.7.8a)
2. collecting an array of alignment metrics using Picard tools (v3.1.1)

### Read mapping - `modules/star/align/main.nf`

>[!NOTE]
>*Besides the STAR mapper algorithm, Samtools is used for processing and analysing the data. It includes tools for file format conversion and manipulation, sorting, querying, and statistics amongst other methods.*

We align the RNA-seq reads to a genome with STAR in order to identify single nucleotide variants (SNVs) and chimeras (fusion genes).
```groovy
process STAR_ALIGN {    
    // Directives
    // input:
    // output:
    script:
    """
    STAR \\
		--runThreadN "$task.cpus" \\
		--genomeDir "${star_index}" \\
		--readFilesCommand "gunzip -c" \\
		--readFilesIn ${trim1} ${trim2} \\
		--outFileNamePrefix "${sample}-" \\
		${star_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR: \$(STAR --version | sed 's/STAR_//')
    END_VERSIONS
    """
}
```
Transcribed to:
```bash        
# e.g. VI-3429-593-2DH-1
STAR \
    --runThreadN 16 \
    --genomeDir ref_genome.star.idx \
    --readFilesCommand "gunzip -c" \
    --twopassMode Basic \
    --readFilesIn VI-3429-593-2DH-1_R1_001_val_1.fq.gz VI-3429-593-2DH-1_R2_002_val_2.fq.gz \
    --outFileNamePrefix VI-3429-593-2DH-1 \
    --chimOutType Junctions \ # **essential** to create the Chimeric.junction.out file for STAR-Fusion
	--chimSegmentMin 12 \
	--chimJunctionOverhangMin 8 \  # **essential to invoke chimeric read detection & reporting**
	--chimOutJunctionFormat 1 \ # **essential** includes required metadata in Chimeric.out.junction file
	--outSAMtype BAM SortedByCoordinate \
	--outBAMcompression 6 \
	--outSAMstrandField intronMotif \
	--outSAMunmapped Within \
	--outReadsUnmapped None \
	--outFilterScoreMinOverLread 0 \
	--outFilterMatchNminOverLread 0 \
	--outFilterMatchNmin 0 \
	--outFilterMismatchNmax 2 \
	--alignSJstitchMismatchNmax 5 -1 5 5 \
	--alignIntronMin 10 \
	--alignIntronMax 100000 \
	--alignMatesGapMax 100000
```
During the alignment step the following **crucial** parameters are applied:
- `--twopassMode Basic` – enambles the most sensitive novel splice junction discovery: during the 1st pass STAR maps with the usual parameters, then in hte 2nd pass it collects the previously detected junctions and uses the mas "annotated" junctions, what allows to detect more spliced reads mapping to novel junctions, 
- `--chimOutType Junctions` – outputs the chimeric alignments to a separate file, called *"Chimeric.out.junction"*, that is the input file for STAR-Fusion,
- `--chimOutJunctionFormat 1` – is essential to include required metadata in Chimeric.junction.out file for STAR-Fusion,
- `--chimSegmentMin [INT]` – detect fusions that map to two chromosomes, with a minimum length above [INT] on either chromosome, this parameter is essential to start chimeric read detection and reporting

> [!NOTE]
> The STAR reference index is precompiled in the CTAT-genome-lib, but as a failsafe we also included the reference generating script `STAR --runMode genomeGenerate` with the same `get_index_channel()` function, we used for Salmon to check missingness. 

### Alignment QC - `modules/picard/rnametrics/main.nf`

Tools, namely the `CollectRnaSeqMetrics` command from the [Picard tools](https://broadinstitute.github.io/picard/) suite (v3.1.1) was used for QC after mapping. 
We produced statistics to gather the necessary information about the success of the alignment, checking the quality of the mapping process: the percentage of mapped reads, the ratio of multi-mapped reads and reads that passed machine signal-to-noise threshold quality filters. These global indicators show the overall sequencing accuracy and of the presence of contaminating DNA. This tool also calculates the total numbers and the fractions of nucleotides within specific genomic regions including untranslated regions (UTRs), introns, intergenic sequences (between discrete genes), and peptide-coding sequences (exons).
```groovy
process PICARD_RNAMETRICS {  
    // Directives
    // input:
    // output:
    script:
    """
    # Collect RNA-seq metrics
    java -jar \$PICARD CollectRnaSeqMetrics \\
        -I ${bam_file} \\
        -O "${sample}.rna_metrics.txt" \\
        --REF_FLAT ${ref_flat} \\
        --REFERENCE_SEQUENCE ${fasta} \\
        --RIBOSOMAL_INTERVALS ${rrna_intervals} \\
        ${params.picard_metrics_args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(java -jar \$PICARD CollectRnaSeqMetrics 2>&1) | grep -oP "(?<=Version:).+?(?=\\ )")
    END_VERSIONS
    """
}
```
Transcribed to:
```bash
# e.g. VI-3429-593-2DH-1
java -jar $PICARD CollectRnaSeqMetrics \
    -I "VI-3429-593-2DH-1-Aligned.SortedByCoord.bam" \
    -O "VI-3429-593-2DH-1.rna_metrics.txt" \
    --REF_FLAT "ref_annot.refFlat.txt" \
    --REFERENCE_SEQUENCE "ref_genome.fa" \
    --RIBOSOMAL_INTERVALS "rRNA_intervals.txt" \
    --STRAND SECOND_READ_TRANSCRIPTION_STRAND \
    --ASSUME_SORTED true
```
> [!IMPORTANT]
> On its input, this tools requires valid SAM/BAM file(s), a tab-delimited file REF_FLAT file containing information about the location of other transcripts, exon start and stop sites, etc., and (optinally) a file including the location of rRNA sequences in the genome, in interval_list format. 

This last two files are not available by default in the CTAT-genome-lib, so we created them, using custom scripts:
```bash
# Sequence names and lengths. (Must be tab-delimited.)
cut -f1,2 "${params.genome_file}.fai" | \
    perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:GRCh38"' | \
    grep -v _ >> "rRNA_intervals.txt"
# Intervals for rRNA transcripts.
grep 'gene_type "rRNA"' "${params.gtf_annot}" | \
    awk '$3 == "transcript"' | \
    cut -f1,4,5,7,9 | \
    perl -lane '/transcript_id "([^"]+)"/ or die "no transcript_id on $.";
                print join "\t", (@F[0,1,2,3], $1)' | \
    sort -k1V -k2n -k3n >> "rRNA_intervals.txt"

# Convert GTF to genePred format
gtfToGenePred \\
    -genePredExt \\
    -geneNameAsName2 \\
    "${params.gtf_annot}" \\
    "ref_annot.refFlat.tmp.txt"
# Parse genePred file to UCSC refFlat format
paste <(cut -f 12 ref_annot.refFlat.tmp.txt) <(cut -f 1-10 ref_annot.refFlat.tmp.txt) > ref_annot.refFlat.txt
rm ref_annot.refFlat.tmp.txt
gzip ref_annot.refFlat.txt
```

## 3.3. Fusion workflow

This subworkflow integrates:
1. Fusion detection using STAR-fusion (v1.10.1)
2. Validation using FusionInspector
3. And *de novo* fusion transcript reconstruction using Trinity (v2.14.0)

## 2.4. Fusion detection - `modules/star_fusion/main.nf`

Following the alignment we used [STAR-fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki) (v1.10.1) to identify candidate fusion transcripts - chimeras - that align to two chromosomes, each segment longer than the threshold set with `--chimSegmentMin`. Since, STAR aligner was run *a priori*, we ran STAR-fusion in "**Kickstart mode**" using the existing  'Chimeric.junction.out' output files. 

> [!NOTE]
> Since we experienced trouble with running STAR-fusion in 'kickstart mode' according to the provider's instructions, I used bith the junction file and the raw read inputs, which solved this issue.

Predicted fusions were then 'in silico validated' using [FusionInspector](https://github.com/FusionInspector/FusionInspector/wiki/installing-FusionInspector), which performs a more refined exploration of the candidate fusion transcripts, runs [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) to *de novo* assemble fusion transcripts from the RNA-Seq reads, and provides the evidence in a suitable format to facilitate visualization. Predicted fusions are annotated according to prior knowledge of fusion transcripts relevant to cancer biology (or previously observed in normal samples and less likely to be relevant to cancer), and assessed for the impact of the predicted fusion event on coding regions, indicating whether the fusion is in-frame or frame-shifted along with combinations of domains that are expected to exist in the chimeric protein.

```groovy
process FUSION {
    // Directives
    // input:
    // output:
    script:
    """
    STAR-Fusion \\
        --genome_lib_dir ${genome_lib_dir} \\
        --left_fq ${trim1} \\
        --right_fq ${trim2} \\
        -J ${junction} \\
        --CPU ${task.cpus} \\
        --output_dir . \\
        ${params.fusion_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """
}
```
Transcribed to:
```bash
# e.g. VI-3429-593-2DH-1
STAR-Fusion \ 
  --genome_lib_dir ctat-genome-lib \ 
  --chimeric_junctions VI-3429-593-2DH-1-Chimeric.junction.out \ 
  --FusionInspector validate \ 
  --denovo_reconstruct \ 
  --examine_coding_effect \ 
  --output-dir VI-3429-593-2DH-1
```
Running parameters:
- `--genome_lib_dir [PATH]` – path to the CTAT genome lib,
- `-J/--chimeric_junctions [PATH]` – accession path to the 'Chimeric.out.junction' file created by STAR
- `--FusionInspector 'validate'` – requires the candidates (in a format: geneA--geneB) from the first column of the 'Chimeric.out.junction' file, then aligns them to the genome to identify those reads that align concordantly as fusion evidence to the fusion contigs,
- `--denovo_reconstruct` – *de novo* reconstruction of fusion transcripts using Trinity, creates a `.fasta` and a `.bed` file,
- `--examine_coding_effect` –  explores the effect the fusion events have on coding regions of the fused genes.

## 3.4. Variant workflow

This subworkflow integrates:
1. Marking PCR duplicated with Picard tools (v3.1.1) using `MarkDuplicates`
2. Subsampling (25%, 50%, 75%) reads at random with Picard tools using `DownsampleSam`
3. Somatic variant calling with BCFtools (v1.20) using `mpileup` and `call`
4. Annotating variants using VEP (v113.0)

While [GATK](https://gatk.broadinstitute.org/hc/en-us) suite is perhaps the most popular variant caller, in my experience it is most sensitive for detecting germline single nucleotide variants (SNVs), and doesn't work as well for RNA-seq. Plus, it requires multiple extra steps of pre-processing. Hence, we worked with the 'old faithful' pipeline, based on [BCFtools](https://samtools.github.io/bcftools/howtos/variant-calling.html), which we found has greater sensitivity for SNVs when compared to GATK. 

> [!NOTE]
> Although, false positive rate can be higher too, and mpileup is neither fine-tuned to identify indels, these limitations we chose to compensate with cross-referencing RNA-seq and WGS results.

#### Remove low MAPQ reads, and mark and remove PCR and/or optical duplicates - `modules/picard/dedup/main.nf`
```groovy
process PICARD_DEDUP {
    // Directives
    // input:
    // output:
    script:
    """
    # Collect RNA-seq metrics
    java -jar \$PICARD MarkDuplicates \\
		--INPUT ${bam_file} \\
		--OUTPUT ${sample}-dedup.bam \\
		--METRICS_FILE ${sample}-dup-metrics.txt \\
		--ASSUME_SORT_ORDER coordinate
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
	picard: \$(echo \$(java -jar \$PICARD MarkDuplicates 2>&1) | grep -oP "(?<=Version:).+?(?=\\ )")
    END_VERSIONS
    """
}
```
#### Produce 3 random BAM subsets per sample (at 25%, 50%, and 75%) - `modules/picard/downsample/main.nf`
```groovy
process PICARD_DOWNSAMPLE {
    // Directives
    // input:
    // output:
    script:
    """
    # Create a 25% subsample of the input BAM file
    java -jar \$PICARD DownsampleSam \\
      --INPUT ${bam_file} \\
      --OUTPUT "${sample}-25pc.bam" \\
      --RANDOM_SEED 42 \\
      --PROBABILITY 0.25 \\
      --VALIDATION_STRINGENCY SILENT

    # Create a 50% subsample of the input BAM file
    java -jar \$PICARD DownsampleSam \\
      --INPUT ${bam_file} \\
      --OUTPUT "${sample}-50pc.bam" \\
      --RANDOM_SEED 42 \\
      --PROBABILITY 0.50 \\
      --VALIDATION_STRINGENCY SILENT

    # Create a 75% subsample of the input BAM file
    java -jar \$PICARD DownsampleSam \\
      --INPUT ${bam_file} \\
      --OUTPUT "${sample}-75pc.bam" \\
      --RANDOM_SEED 42 \\
      --PROBABILITY 0.75 \\
      --VALIDATION_STRINGENCY SILENT
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(java -jar \$PICARD DownsampleSam 2>&1) | grep -oP "(?<=Version:).+?(?=\\ )")
    END_VERSIONS
    """
}
```
#### Call variants - `modules/bcftools/mpileup/main.nf`

We called variants using BCFtools `mpileup` on the main BAM and 3 subsets together, piped into BCFtools `call`. The first mpileup part generates genotype likelihoods at each genomic position with coverage. The second call part makes the actual calls. 

```groovy
process BCFTOOLS_MPILEUP {
    // Directives
    // input:
    // output:
    script:
    """
    bcftools mpileup \\
        --threads "${task.cpus}" \\
        -Ou --fasta-ref "${ref_genome}" \\
        ${bam} ${bam_sub1} ${bam_sub2} ${bam_sub3} \\
        ${params.snv_args_mpileup} |\\
        bcftools call \\
        ${params.snv_args_call} > "${sample}.vcf.gz"

    tabix -p vcf -f ${sample}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) |  grep -oP "(?<=Version: ).+?(?=\\ )")	
    END_VERSIONS
    """
}
```
Transcribed to:
```bash
# e.g. VI-3429-593-2DH-1
bcftools mpileup \\
    --threads 16 \\
    -Ou --fasta-ref ref_genome.fa \\
    VI-3429-593-2DH-1-dedup.bam VI-3429-593-2DH-1-25pc.bam VI-3429-593-2DH-1-50pc.bam VI-3429-593-2DH-1-75pc.bam \\
    --redo-BAQ \\
    --min-BQ 30 \\
    --per-sample-mF \\
    --annotate "FORMAT/DP,FORMAT/AD" |\\
    bcftools call \\
    --multiallelic-caller \\
    --variants-only > VI-3429-593-2DH-1.vcf.gz
```
Running parameters `bcftools mpileup`:
- `-Ou` – outputs uncompressed BCF format to standard output (u for uncompressed, O for BCF), ideal for piping into bcftools call.
- `--redo-BAQ` – recomputes Base Alignment Quality (BAQ), which helps reduce false positives near indels
- `--min-BQ 30` – filters out bases with base quality scores below 30, ensuring high-confidence variant calls.
- `--per-sample-mF` – enables per-sample filtering of reads flagged as secondary, supplementary, or failing vendor quality checks: `m` corresponds to the min number of gapped reads for indel candidates (default: 1) and `F` corresponds to the min fraction of gapped reads (default: 0.002)
- `--annotate "FORMAT/DP,FORMAT/AD"` – adds depth of coverage (DP) and allele depth (AD) annotations to the VCF output, providing insight into allele support per sample

Running parameters `bcftools call`:
- `--multiallelic-caller` – enables multiallelic calling model, allowing detection of variants with more than two alleles at the same position.
- `--variants-only` – excludes reference sites from the output, reporting only positions with observed variants.

# 4. Complete workflow

Is the entry workflow, which any script can define up to one of. The complete workflow does not have a name and serves as the entrypoint of the script, that we call from the command line.

```groovy
#!/usr/bin/env nextflow
/*
 * Import processes from modules
 */
include { core_workflow } from './subworkflows/core-workflow.nf'
include { alignment_workflow } from './subworkflows/alignment-workflow.nf'
include { fusion_workflow } from './subworkflows/fusion-workflow.nf'
include { variant_workflow } from './subworkflows/variant-workflow.nf'
include { MULTIQC } from './modules/multiqc/main.nf'

process COMPILE_VERSIONS {
    input:
    path version_files

    output:
    path "versions.yml"

    script:
    """
    echo "---" > versions.yml
    for f in ${version_files}; do
    	cat \$f >> versions.yml
    done
    """
}

workflow {
    /*
     * Instantiate empty output channels
     */
    report_ch = Channel.empty()
    version_ch = Channel.empty()
    /*
     * RNA-Seq Core Workflow
     */
    // Read in read fastq files
    read_pairs_ch = Channel.fromPath(params.reads, checkIfExists: true) \
        | splitCsv(header: true, sep: ',', strip: true) \
        | map { row -> tuple(row.sample, file(row.read1), file(row.read2)) }
    // Run core workflow
    core_workflow(read_pairs_ch)
    trimmed_reads = core_workflow.out.trimmed_reads
    report_ch.mix(core_workflow.out.report)
    /*
     * RNA-Seq Alignment Workflow
     */
    alignment_workflow(trimmed_reads)   
    bam_files = alignment_workflow.out.bam
    junctions = alignment_workflow.out.junction
    report_ch.mix(alignment_workflow.out.report)
    /*
     * Gene-fusion Detection Workflow
     */
    fusion_workflow(trimmed_reads, junctions)
    /*
     * Single Nucleotid Variant Calling Workflow
     */
    variant_workflow(bam_files)
    report_ch.mix(variant_workflow.out.report)
    /*
     * MultiQC
     */
    MULTIQC(report_ch.collect())
    /*
     * Save module versions
     */
    version_ch.mix(
	    core_workflow.out.version,
	    alignment_workflow.out.version,
	    //fusion_workflow.out.version,
	    variant_workflow.out.version).collect()

    COMPILE_VERSIONS(version_ch)
    
}
```
In the main workflow the `Channel.fromPath()` channel factory creates the input channel from the `samples.csv` file. The file rows are read in as strings, separated using `splitCsv(header: true, sep: ',', strip: true)` channel operator, and formatted into the required tuple using `map { row -> tuple(row.sample, file(row.read1), file(row.read2)) }`. 

#### Generating a report - `modules/multiqc/main.nf`

As a last step we produce a report by gathering and parsing summary statistics from results and log files generated by all other bioinformatics tools with **MultiQC (v.1.10.1)**. MultiQC recursively searches through any provided file paths and finds files that it recognises. 

It parses relevant information from these and generates a single stand-alone HTML report file. Furthermore, a multiqc_data folder will be generated as well with the reformatted compact data tables in there, one from each module in TSV format, as well as a verbose multiqc.log file and a multiqc_data.json.

Whilst MultiQC is typically used as a final reporting step in an analysis, it can also be used as an intermediate in your analysis, as these files essentially standardize the outputs from a lot of different tools and make them useful for downstream analysis. We will do that as well, going one step further with the TidyMultiqc package. The TidyMultiqc package provides the means to convert the multiqc_data.json file into a tidy data frame for downstream analysis in R. 

## 4.1. Running the Pipeline
```bash
nextflow run complete-workflow.nf -C nextflow.config -profile uppmax --reads $NXF_HOME/meta/samples.csv --reference $NXF_HOME/ctat-genome-lib -with-report $NXF_HOME/report.html -with-timeline $NXF_HOME/timeline.html
```

With these specification Nextflow will create an HTML execution report and another HTML timeline. The execution report includes three main sections: `summary`, `resources` and `tasks`. The summary section reports the execution status, the launch command, overall execution time and some other workflow metadata. The resources section plots the distribution of resource usage for each workflow process using the interactive plotly.js plotting library. Finally, the tasks section lists all executed tasks, reporting for each of them the status, the actual command script, and many other metrics. The rendered timeline includes all executed processes from our pipeline, where each bar represents a process run in the pipeline execution with additional textual information about the task duration time and the virtual memory size peak.

> [!IMPORTANT]
>***When running the nextflow pipeline it is IMPORTANT (!!!) to also move the working directory to $NXF_HOME. Because nextflow automatically outputs the `work/` directory to `pwd` and it can overload the allocated resource on the entry node.***