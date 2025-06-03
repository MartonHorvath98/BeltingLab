process FASTQC {
    tag "FastQC on $sample"
    label "low_effort"

    publishDir "${params.report_outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)
    val subfolder

    output:
    tuple val(sample), path("$subfolder/${sample}/*_fastqc.zip"), emit: fastqc_report
    path "versions.yml"                                            , emit: version

    script:
    """
    mkdir -p $subfolder/${sample}
    fastqc \\
        -o $subfolder/${sample} \\
        ${params.fastqc_args} \\
        --threads ${task.cpus} \\
        ${read1} ${read2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """

}
