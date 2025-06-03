#!/usr/bin/env nextflow

/*
 * Import processes from modules
 */
include { FUSION } from '../modules/star_fusion/main.nf'


workflow fusion_workflow {
    /*
     * Check if the reference genome and STAR index are provided
     */
    take:
    trimmed_reads
    junctions

    main:
    log.info """\
    R N A S E Q - N F   F U S I O N - D E T E C T I O N
    =====================================================
    ctat-genome-lib reference : ${params.reference}
    STAR alignment output     : ${params.star_outdir}
    =====================================================
    """
    .stripIndent(true)
    /*
     * Instantiate empty output channels
     */
    version_ch = Channel.empty()
    /*
     * STAR-Fusion 
     */
    // Load CTAT reference 
    genome_lib_ch = Channel.value(params.reference)
    // Create input channel
    input_ch = trimmed_reads.join(junctions).map { tuple ->
        def sample_id = tuple[0]
        def trim1 = tuple[1]
        def trim2 = tuple[2]
        def junction = tuple[3]
        return [sample_id, trim1, trim2, junction]
    }
    // Run fusion detection
    FUSION(genome_lib_ch, input_ch)
    /*
     * Collect version informations
     */
    version_ch.mix(FUSION.out.version).collect()

    emit:
    fusion = FUSION.out.fusions // channel: [ val(sample), path(fusion.tsv) ] 
    version = version_ch
}

