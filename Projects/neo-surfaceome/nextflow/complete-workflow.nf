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
    //fusion_workflow(trimmed_reads, junctions)
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

