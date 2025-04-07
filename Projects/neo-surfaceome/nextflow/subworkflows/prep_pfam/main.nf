#!/usr/bin/env nextflow
/*
 * Nextflow module to download and extract Pfam database files.
 * The module uses wget to download files and stores them in the specified output directory.
 * The module uses gunzip to extract the downloaded files.
 */

// Pipeline parameters
params.url="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/"
params.filenames=["Pfam-A.hmm.gz", "active_site.dat.gz"]
params.dataset="pfam"

// Include modules
include { GET } from '../../modules/get/main.nf'
include { HMMPRESS } from '../../modules/hmmer/hmmpress/main.nf'

workflow {
    // Create the channel with filenames
    file_ch = Channel.of(params.filenames).flatten()
    // STEP 1: Download the files
    download_ch = GET(file_ch)
    // STEP 2: Select out the HMM file
    hmm_file_ch = download_ch.filter { it.name.endsWith('.hmm') }
    hmm_file_ch.view { file ->
        "HMM file: $file"
    }
    // STEP 3: Run HMMER HMMpress on the HMM file
    // This will create the .h3m, .h3i, .h3f, and .h3p files
    HMMPRESS(hmm_file_ch)
    
}
