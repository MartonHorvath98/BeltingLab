#!/usr/bin/env nextflow
/*
* Nextflow module to dowload files from //https:, //ftp:, //http: host sites on the internet.
* The module uses wget to download files and stores them in the specified output directory, then unzips them using gunzip.
*/
process GET {
    publishDir "${params.ref_dir}", mode: 'copy'
    
    input:
    val url

    output:
    path input_file, emit: archive

    script:
    input_file = url.split('/').last()
    """
    wget -O $input_file $url
    """
}