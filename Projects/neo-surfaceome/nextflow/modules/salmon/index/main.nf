#!/usr/bin/env nextflow
/*
 * Nextflow module to create a salmon index from a genome and transcriptome fasta file.
 * The script creates a gentrome, by concatenating the genome and transcriptome fasta files, 
 * for more accurate quantification of the transcriptomic data. 
*/
process SALMON_INDEX{
    publishDir "${params.base_publish_dir}/${params.dataset}", mode: 'copy'
    label params.dataset
    
    input:
    tuple path(genome_fasta), path(transcriptome_fasta)

    output:
    path index_folder, emit: salmon_index

    script:
    index_folder = file("${genome_fasta.baseName}salmon.idx")
    """
    # Create decoy file for salmon index
    grep "^>" ${genome_fasta} | cut -d " " -f 1 | sed -e 's/>//g' > "${params.base_publish_dir}/${params.dataset}/decoys.txt"
    # Concatenate transcriptome and genome
    cat ${transcriptome_fasta} ${genome_fasta} > "${params.base_publish_dir}/${params.dataset}/gentrome.fa"
    # Create salmon index
    salmon index -t ${params.base_publish_dir}/${params.dataset}/gentrome.fa -d ${params.base_publish_dir}/${params.dataset}/decoys.txt \
    -i ${index_folder} --gencode -k 31 --threads ${task.cpus}
    # Remove intermediate files
    rm "${params.base_publish_dir}/${params.dataset}/gentrome.fa"
    rm "${params.base_publish_dir}/${params.dataset}/decoys.txt"
    """
}