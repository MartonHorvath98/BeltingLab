#!/usr/bin/env nextflow

/*
 * Import processes from modules
 */
include { TRIM_GALORE } from './modules/trim_galore/main.nf'
// include { MULTIQC} from './modules/multiqc/main.nf'
include { INDEX } from './modules/salmon/index/main.nf'
// include { QUANTIFICATION } from './modules/salmon/quant/main.nf'
/*
 * pipeline input parameters
 */

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
    trim_ch = TRIM_GALORE(read_pairs_ch)
    
    /*
     * QUANTIFICATION USING SALMON
     */
    // Build index for salmon
    index_ch = INDEX(params.transcriptome_file, params.genome_file)
    // fastqc_ch = FASTQC(read_pairs_ch)
    // Create MultiQC report
    // MULTIQC(quant_ch.mix(fastqc_ch).collect())

}

