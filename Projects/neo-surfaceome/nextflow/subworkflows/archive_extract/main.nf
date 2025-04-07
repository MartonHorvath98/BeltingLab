#!/usr/bin/env nextflow
/*
 * Nextflow module to download and extract CTAT database files.
 */

// Include modules
include { UNTAR } from '../../modules/untar/main.nf'
include { GUNZIP } from '../../modules/gunzip/main.nf'

// Custom process to move the files to the output directory
// Process to extract the CTAT database files
workflow ARCHIVE_EXTRACT {
    take:
    archive
 
    main:
    archive.view { "Archive input: ${it}" }
    
    archive
    .branch {
        tar: it.toString().endsWith('.tar.gz')
        gz: it.toString().endsWith('.gz')
        non_assigned: true
    }
    .set { archive_to_extract }  

    // This is a confidence check
    not_extracted = archive_to_extract.non_assigned
    not_extracted.view { archive_ -> log.warn("Archive not in an expected format: " + archive_) }
    

    GUNZIP(archive_to_extract.gz)
    UNTAR(archive_to_extract.tar)

    extracted = Channel
            .empty()
            .mix(
                GUNZIP.out.dataset,
                UNTAR.out.dataset
            )

    extracted.view { archive_ -> log.info("Extracted archive: " + archive_) }

    emit:
    not_extracted
    extracted

}