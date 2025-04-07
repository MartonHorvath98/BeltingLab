#!/usr/bin/env nextflow
/*
* Nextflow module to extract tar compressed files into a specified output directory.
*/
process UNTAR {
    publishDir "${params.ctat_lib}", mode: 'copy'
    
    input:
    path archive

    output:
    path "${prefix}", emit: dataset

    script:
    prefix = archive.baseName.toString().replaceFirst(/\\.*/,'')
    """
    mkdir -p ${prefix}
    
    ## If just files or multiple directories, place all in prefix
    if [[ \$(tar -taf ${archive} | grep -o -P "^.*?\\/" | uniq | wc -l) -eq 1 ]]; then
        tar \\
            -C ${prefix} --strip-components 1 \\
            -xavf ${archive} \\
            ${args} 
    else
        tar \\
            -C ${prefix} \\
            -xavf ${archive} \\
            ${args}
    fi

    rm -f ${archive}
    """
}