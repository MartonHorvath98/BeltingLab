#!/usr/bin/env nextflow

/*
 * Import processes from modules
 */
include { INDEX } from './modules/star/index/main.nf'
include { ALIGNMENT } from './modules/star/align/main.nf'
include { FUSION } from './modules/star-fusion/main.nf'
/*
 * pipeline functions
 */
def get_index_channel() {
    def index_path = file(params.star_index)
    return index_path.exists()
        ? Channel.value(index_path)
        : INDEX(params.genome_file, params.gtf_annot)
}

workflow {
    log.info """\
    R N A S E Q - N F   A D V A N C E D   P I P E L I N E
    =====================================================
    ctat-genome-lib reference : ${params.reference}
    read meta file            : ${params.reads}
    =====================================================
    STEP 1: ALIGNMENT (STAR)
    STEP 2: FUSION EVENT DETECTION (STAR-FUSION)
    """
    .stripIndent(true)
    // Read in trimmed fastq file pairs from the trim-galore output
    input_ch = Channel.fromFilePairs('results/01_trim/*_val_{1,2}.fq.gz', checkIfExists: true, flatten: true)
                      .map { tuple_row -> 
                            def (sample_id, files) = tuple_row
                            def (read1, read2) = files
                            return tuple(sample_id, read1, read2)
                      }
    /*
     * ALIGNMENT USING STAR
     */
    index_ch = get_index_channel()
    // Align reads using STAR
    align_input_ch = index_ch.combine(input_ch)
                             .view{ "Alignment input: $it" }
    ALIGNMENT(align_input_ch)
    /*
     * FUSION EVENT DETECTION USING STAR-FUSION
     */
    fusion_input_ch = input_ch
    .join(ALIGNMENT.out.junction, by: 0)
    .map { sample_id, trim1, trim2, junction ->
        tuple(params.reference, sample_id, trim1, trim2, junction)
    }
    FUSION(fusion_input_ch)

}

