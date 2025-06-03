#!/usr/bin/env nextflow

/*
 * Import processes from modules
 */
include { FASTQC as PRE_FASTQC } from '../modules/fastqc/main.nf'
include { FASTQC as POST_FASTQC } from '../modules/fastqc/main.nf'
include { SEQKIT as PRE_SEQKIT } from '../modules/seqkit/main.nf'
include { SEQKIT as POST_SEQKIT } from '../modules/seqkit/main.nf'
include { TRIM_GALORE } from '../modules/trim_galore/main.nf'
include { SALMON_INDEX as INDEX } from '../modules/salmon/index/main.nf'
include { SALMON_QUANT as QUANTIFICATION } from '../modules/salmon/quant/main.nf'
/*
 * pipeline functions
 */
def get_index_channel() {
    def index_path = file(params.salmon_index)
    return index_path.exists()
        ? Channel.value(index_path)
        : INDEX(params.transcriptome_file, params.genome_file)
}

workflow core_workflow {
    take:
    read_pairs_ch

    main:
    log.info """\
    R N A S E Q - N F   C O R E - P I P E L I N E
    =============================================
    ctat-genome-lib reference : ${params.reference}
    read meta file            : ${params.reads}
    =============================================
    STEP 1: TRIMMING (TRIM-GALORE)
    STEP 2: QUANTIFICATION (SALMON)
    STEP 3: REPORTING (MULTIQC)
    """
    .stripIndent(true)
    /*
     * Instantiate empty output channels
     */
    version_ch = Channel.empty()
    report_ch = Channel.empty()
    /*
     * RAW READ QUALITY CHECK USING FASTQC AND SEQKIT
     */
    def subfolder = "raw_data"
    PRE_FASTQC(read_pairs_ch, subfolder) // Run fastqc on raw reads
    PRE_SEQKIT(read_pairs_ch, subfolder)
    /*
     * TRIMMING USING TRIM-GALORE
     */
    // Trim low quality and adapter contaminated reads
    TRIM_GALORE(read_pairs_ch)
    /*
     * POST-TRIMMING QUALITY CHECK USING FASTQC AND SEQKIT
     */
    subfolder = "trimmed_data"
    POST_FASTQC(TRIM_GALORE.out.trimmed_reads, subfolder)
    POST_SEQKIT(TRIM_GALORE.out.trimmed_reads, subfolder) // Run seqkit to check for read quality
    /*
     * QUANTIFICATION USING SALMON
     */
    // Build index for salmon
    index_ch = get_index_channel()
    // Quantify reads using salmon
    quant_input_ch = index_ch.combine(TRIM_GALORE.out.trimmed_reads)
    QUANTIFICATION(quant_input_ch)
    /*
     * Collect logs for MultiQC report
     */
    report_ch.mix(TRIM_GALORE.out.log.map{ it[1] },POST_FASTQC.out.fastqc_report.map{ it[1] },QUANTIFICATION.out.quant_files.map{ it[1] }).collect()
    /*
     * Collect version informations
     */
    version_ch.mix(POST_FASTQC.out.version,POST_SEQKIT.out.version,QUANTIFICATION.out.version).collect()

    emit:
    trimmed_reads = TRIM_GALORE.out.trimmed_reads // channel: [ val(sample), path(trim1), path(trim2) ]
    version = version_ch
    report = report_ch
}

