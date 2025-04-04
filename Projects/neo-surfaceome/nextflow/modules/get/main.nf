#!/usr/bin/env nextflow
/*
* Nextflow module to dowload files from //https:, //ftp:, //http: host sites on the internet.
* The module uses wget to download files and stores them in the specified output directory, then unzips them using gunzip.
*/
process GET{
    publishDir "${params.base_publish_dir}/${params.dataset}", mode: 'copy'
    label params.dataset
    
    input:
    val url
    val input_file

    output:
    path output_file

    script:
    output_file = input_file - '.gz'
    """
    wget -O $input_file $url/$input_file 
    gunzip -c $input_file > $output_file
    """
}