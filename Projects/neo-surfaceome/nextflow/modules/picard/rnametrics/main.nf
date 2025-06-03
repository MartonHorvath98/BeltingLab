 process PICARD_RNAMETRICS {
    label "generic_task"
    tag "Collect RNA-seq metrics from ${sample}"
    publishDir "${params.star_outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(bam_file)
    path ref_flat
    path fasta
    path rrna_intervals

    output:
    tuple val(sample), path("*.rna_metrics.txt"), emit: metrics
    tuple val(sample), path("*.pdf")            , emit: pdf, optional: true
    path "versions.yml"                            , emit: version

    script:
    """
    # Collect RNA-seq metrics
    java -jar \$PICARD CollectRnaSeqMetrics \\
        -I ${bam_file} \\
        -O "${sample}.rna_metrics.txt" \\
        --REF_FLAT ${ref_flat} \\
        --REFERENCE_SEQUENCE ${fasta} \\
        --RIBOSOMAL_INTERVALS ${rrna_intervals} \\
        ${params.picard_metrics_args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(java -jar \$PICARD CollectRnaSeqMetrics 2>&1) | grep -oP "(?<=Version:).+?(?=\\ )")
    END_VERSIONS
    """
}
