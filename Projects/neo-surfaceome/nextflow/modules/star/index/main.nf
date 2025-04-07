#!/usr/bin/env nextflow
/* 
 * Nextflow module to create a STAR index from a genome fasta file.
 * The script uses the STAR aligner to generate the index files required for RNA-seq alignment.
 */
process STAR_INDEX {
    publishDir "${params.base_publish_dir}/${params.genome_file}.star.idx", mode: 'copy'
    label params.dataset
    
    conda "${moduleDir}/environment.yml"

    input:
    path genome_fasta
    path gtf_annotation

    output:
    path "star", emit: star_index

    script:
    """
    # Create STAR index
    STAR --runThreadN $params.cpus --runMode genomeGenerate --genomeDir star \\
         --genomeFastaFiles $params.genome_file --sjdbGTFfile $params.gtf_file --sjdbOverhang 100
    """
}