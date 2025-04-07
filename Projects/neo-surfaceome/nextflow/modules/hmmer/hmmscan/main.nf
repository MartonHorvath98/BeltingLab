#!/usr/bin/env nextflow
/*
 * Nextflow module to run hmmscan on a set of protein sequences against a HMM database.
 * The script uses the HMMER tool to perform the scan and outputs the results in a tabular format.
 */
process HMMSCAN{
    publishDir "${params.base_publish_dir}/${params.dataset}", mode: 'copy'
    label params.dataset
    
    input:
    path protein_fasta
    path hmm_database

    output:
    path hmmscan_output, emit: hmmscan_results

    script:
    hmmscan_output = file("${protein_fasta.baseName}.hmmscan.out")
    """
    # Run hmmscan on the protein sequences against the HMM database
    hmmscan ${params.hmmscan_args} --cpu ${task.cpus} --domtblout ${hmmscan_output} ${hmm_database} ${protein_fasta}
    """
}