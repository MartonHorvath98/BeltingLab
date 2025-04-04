#!/usr/bin/env nextflow
/*
 * Nextflow module to download and extract 3DID database files. The module uses wget to 
 * download files, gunzip to extract them and stores them in the specified output directory.  
 * The module also processes the downloaded files to extract interacting domains using a
 * custom perl script.
 */

// Pipeline parameters
params.url="https://3did.irbbarcelona.org/download/current/"
params.filenames=["3did_flat.gz"]
params.dataset="did"

// Include modules
include { GET } from '../../modules/get/main.nf'

process PREPROCESS {
    publishDir "${params.base_publish_dir}/${params.dataset}", mode: 'copy'
    label params.dataset

    input:
    path input_file
    val processed_name

    output:
    path processed_name

    script:
    """
    less "$input_file" | grep "^#=ID" | cut -f4,5 | perl -ane '
    \$F[0] =~ s/.*(PF\\d+).*/\$1/;
    \$F[1] =~ s/.*(PF\\d+).*/\$1/;
    print "\$F[0]\\t\$F[1]\\n\$F[1]\\t\$F[0]\\n";' | sort -u > $processed_name
    """ 
}

workflow {
    // Create the channel with filenames
    file_ch = Channel.of(params.filenames).flatten()
    // STEP 1: Download the files
    download_ch = GET(params.url, file_ch)
    // STEP 2: Filter out the .dat file for processing (uncompressed)
    output_name = download_ch
        .filter { it.name.endsWith('flat') }
        .map { it.simpleName.split('\\_')[0] + '_interactions.tsv' }
        
    // STEP 3: Process domain interactions
    PREPROCESS(download_ch, output_name)
    
}
