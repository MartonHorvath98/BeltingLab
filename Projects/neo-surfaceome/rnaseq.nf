#!/usr/bin/env nextflow

/*
 * Import processes from modules
 */
include { TRIM_GALORE } from './modules/trim_galore/main.nf'
include { MULTIQC} from './modules/multiqc/main.nf'
include { INDEX } from './modules/salmon/index/main.nf'
include { QUANTIFICATION } from './modules/salmon/quant/main.nf'
/*
 * pipeline functions
 */
def get_index_channel() {
    def index_path = file(params.salmon_index)
    return index_path.exists()
        ? Channel.value(index_path)
        : INDEX(params.transcriptome_file, params.genome_file)
}

workflow {
    log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    genome       : ${params.genome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent(true)
    // Read in read fastq files
    read_pairs_ch = Channel.fromPath(params.reads) \
        | splitCsv(header: true, sep: ',', strip: true) \
        | map { row -> tuple(row.sample, file(row.read1), file(row.read2)) }
    /*
     * TRIMMING USING TRIM-GALORE
     */
    // Trim low quality and adapter contaminated reads
    TRIM_GALORE(read_pairs_ch)
    /*
     * QUANTIFICATION USING SALMON
     */
    // Build index for salmon
    index_ch = get_index_channel()
    // Quantify reads using salmon
    quant_input_ch = index_ch.combine(TRIM_GALORE.out.trimmed_reads)
    QUANTIFICATION(quant_input_ch)
    // Create MultiQC report
    report_ch = TRIM_GALORE.out.fastqc_report.mix(QUANTIFICATION.out.quant_files.map{ it[1] })
    report_ch.view { it ->
        log.info "MultiQC report will be generated from ${it.size()} files"
    }.collect().view{
        log.info "MultiQC report will be generated from: ${it}"
    }
    MULTIQC(report_ch.collect())
}

