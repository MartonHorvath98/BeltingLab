#!/usr/bin/env nextflow

/*
 * Import processes from modules
 */
include { PICARD_DEDUP as DEDUP } from '../modules/picard/dedup/main.nf'
include { PICARD_DOWNSAMPLE as DOWNSAMPLE } from '../modules/picard/downsample/main.nf'
include { BCFTOOLS_MPILEUP as MPILEUP } from '../modules/bcftools/mpileup/main.nf'
include { BCFTOOLS_STATS as STATS } from '../modules/bcftools/stats/main.nf'
include { VEP_ANNOTATE as ANNOTATE } from '../modules/vep/main.nf'

workflow variant_workflow {
    /*
     * Check if the reference genome and STAR index are provided
     */
    take:
    bam_files

    main:
    log.info """\
    R N A S E Q - N F   V A R I A N T - C A L L I N G
    =====================================================
    ctat-genome-lib reference : ${params.reference}
    STAR alignment output     : ${params.star_outdir}
    =====================================================
    STEP 1: DEDUPLICATE BAM FILE (PICARD)
    STEP 2: SUBSAMPLE BAM FILE (PICARD)
    STEP 3: MPILEUP & CALL (BCFTOOLS)
    STEP 4: COLLECT STATS (BCFTOOLS)
    STEP 5: ANNOTATE VARIANTS (VEP)
    STEP 6: REPORTING (MULTIQC)
    =====================================================
    """
    .stripIndent(true)
    /*
     * Instantiate empty output channels
     */
    variant_ch = Channel.empty()
    report_ch = Channel.empty()
    version_ch = Channel.empty()
    /*
     * Variant calling (BCFtools)
     */
    // Load reference genome fasta file
    genome_file = Channel.value(params.genome_file)
    // DEDUPLICATE BAM FILE (PICARD)
    DEDUP(bam_files)
    DOWNSAMPLE(DEDUP.out.bam)
    // Create input channel
    input_ch = DEDUP.out.bam.join(DOWNSAMPLE.out.sub_bams).map { tuple ->
        def sample_id = tuple[0]
        def bam = tuple[1]
        def bam_sub1 = tuple[2]
        def bam_sub2 = tuple[3]
        def bam_sub3 = tuple[4]
        return [sample_id, bam, bam_sub1, bam_sub2, bam_sub3]
    }
    // Generating genotype likelyhoods & calling variants
    MPILEUP(genome_file, input_ch)
    /*
     * COLLECT STATS (BCFTOOLS)
     */
    STATS(genome_file, MPILEUP.out.vcf)
    /*
     * ANNOTATE VARIANTS (VEP)
     */
    ANNOTATE(genome_file, MPILEUP.out.vcf)
    // Collect outputs
    /*
     * Collect logs for MultiQC report
     */
    report_ch.mix(DEDUP.out.metrics.map{ it[1] }, ANNOTATE.out.report.map{ it[1] }, STATS.out.stats.map{ it[1] }).collect()
    /*
     * Collect version informations
     */
    version_ch.mix(DEDUP.out.version, DOWNSAMPLE.out.version, MPILEUP.out.version, STATS.out.version, ANNOTATE.out.version).collect()

    emit:
    variant = MPILEUP.out.vcf
    annot = ANNOTATE.out.tab
    report = report_ch
    version = version_ch 
}

