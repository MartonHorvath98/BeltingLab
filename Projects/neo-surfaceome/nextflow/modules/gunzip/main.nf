#!/usr/bin/env nextflow
/*
* Nextflow module to extract gzipped files into a specified output directory.
*/
process GUNZIP {
    publishDir "${params.ctat_lib}", mode: 'copy'
    
    input:
    path archive

    output:
    path "${prefix}", emit: dataset

    script:
    prefix = archive.baseName.toString() - '.gz'
    """
    gunzip  ${prefix} ${archive}
    """
}