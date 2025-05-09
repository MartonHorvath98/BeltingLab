#!/usr/bin/env nextflow

/*
 * Import processes from modules
 */
include { EXTRACT_RRNA } from './modules/other/extract_rrna/main.nf'
include { AGAT_CONVERT } from './modules/other/agat_gtf2bed/main.nf'

workflow {
    log.info """\
    R N A S E Q - N F   A L I G N M E N T   Q C 
    =====================================================
    ctat-genome-lib reference : ${params.reference}
    alignment files (.bam)    : ${params.star_outdir}
    =====================================================
    STEP 1: CONVERT GTF-2-BED (AGAT)
    STEP 2: EXTRACT rRNA SEQUENCES
    STEP 3: COLLECT RNA-SEQ METRICS (PICARD)
    """
    .stripIndent(true)
    // Read in read bam files from folder
    bam_files_ch = Channel.fromPath("${params.star_outdir}/*.bam", checkIfExists: true)
        | map { file(it) }
        | map { file -> tuple(file.name, file) }

    bam_files_ch.view { "BAM files: $it" }

    // Convert GTF to BED format
    bed_annot_ch = AGAT_CONVERT()
    // Make an interval_list file suitable for CollectRnaSeqMetrics.jar
    rrna_intervals_ch = EXTRACT_RRNA()
    rrna_intervals_ch.view { "rRNA intervals: $it" }


    // Align reads using STAR
    // align_input_ch = index_ch.combine(input_ch)
    //                         .view{ "Alignment input: $it" }
    // ALIGNMENT(align_input_ch)
    /*
     * FUSION EVENT DETECTION USING STAR-FUSION
     */
    // fusion_input_ch = input_ch
    // .join(ALIGNMENT.out.junction, by: 0)
    // .map { sample_id, trim1, trim2, junction ->
    //     tuple(params.reference, sample_id, trim1, trim2, junction)
    // }
   // FUSION(fusion_input_ch)

}

