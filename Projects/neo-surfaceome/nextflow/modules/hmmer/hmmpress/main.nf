#!/usr/bin/env nextflow
/*
 * Nextflow module to prepare an HMM database for faster hmmscan searches 
 * using the hmmpress command from the HMMER package (version 3.3.2, Nov 2020)
 */

process HMMPRESS {
    publishDir "${params.base_publish_dir}/${params.dataset}", mode: 'copy'
    label params.dataset

    conda "${moduleDir}/environment.yml"

    input:
    path hmm_file

    output:
    path "${hmm_file}.h3m"
    path "${hmm_file}.h3i"
    path "${hmm_file}.h3f"
    path "${hmm_file}.h3p"

    script:
    """
    hmmpress $hmm_file
    """
}
