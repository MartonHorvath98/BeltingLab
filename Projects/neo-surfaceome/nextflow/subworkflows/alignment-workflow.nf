#!/usr/bin/env nextflow

/*
 * Import processes from modules
 */
include { STAR_INDEX as INDEX } from '../modules/star/index/main.nf'
include { STAR_ALIGN as ALIGNMENT } from '../modules/star/align/main.nf'
include { EXTRACT_RRNA } from '../modules/other/extract_rrna/main.nf'
include { GTF2REFFLAT } from '../modules/other/ucsc_refflat/main.nf'
include { PICARD_RNAMETRICS as RNAMETRICS } from '../modules/picard/rnametrics/main.nf'

/*
 * pipeline functions
 */
def get_index_channel() {
    def index_path = file(params.star_index)
    return index_path.exists()
        ? Channel.value(index_path)
        : INDEX(params.genome_file, params.gtf_annot)
}

workflow alignment_workflow {
    /*
     * Check if the reference genome and STAR index are provided
     */
    take:
    trimmed_reads

    main:
    log.info """\
    R N A S E Q - N F   A L I G N M E N T   Q C 
    =====================================================
    ctat-genome-lib reference : ${params.reference}
    STAR index               : ${params.star_index}
    alignment output (.bam)    : ${params.star_outdir}
    =====================================================
    STEP 1: ALIGN TRIMMED READS (STAR)
    STEP 2: CREATE REFFLAT FILE (UCSC)
    STEP 3: COLLECT RNA-SEQ METRICS (PICARD)
    STEP 4: REPORTING (MULTIQC)
    =====================================================
    """
    .stripIndent(true)
    /*
     * Instantiate empty output channels
     */
    version_ch=Channel.empty()
    report_ch = Channel.empty()
    /*
     * STAR alignment
     */
    // Create STAR index if missing
    index_ch = get_index_channel()
    // STAR alignment
    align_input_ch = index_ch.combine(trimmed_reads)
    ALIGNMENT(align_input_ch)
    /*
     * Extract rRNA sequences from the genome
     */
    EXTRACT_RRNA()
    /*
     * Convert GTF file to REFFlat format
     */
    GTF2REFFLAT()
    /*
     * Collect RNA-seq metrics
     */
    RNAMETRICS(ALIGNMENT.out.out_bam, GTF2REFFLAT.out.ref_flat, params.genome_file, EXTRACT_RRNA.out.rrna_intervals)
    /*
     * Collect logs for MultiQC report
     */
    report_ch.mix(RNAMETRICS.out.metrics.map{ it[1] }, ALIGNMENT.out.log.map{ it[1] }).collect()
    /*
     * Collect version informations
     */
    version_ch.mix(
        ALIGNMENT.out.version,
	EXTRACT_RRNA.out.version,
	GTF2REFFLAT.out.version,
	RNAMETRICS.out.version
	).collect()	

    emit:
    bam = ALIGNMENT.out.out_bam // channel: [ val(sample), path(bam) ]
    junction = ALIGNMENT.out.junction // channel: [ val(sample), path(out.junction) ]
    report = report_ch
    version = version_ch
}

