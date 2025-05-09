process FASTQC {
    tag "FastQC on $sample_id"
    label "low_effort"

    publishDir "${params.report_outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(read1), path(read2)
    val subfolder

    output:
    path "versions.yml"                         , emit: versions
    path "$subfolder/${sample_id}/*_fastqc.zip" , emit: fastqc_report

    script:
    def memory = "${task.memory.toMega()}"
    """
    mkdir -p $subfolder/${sample_id}
    fastqc \\
        -o $subfolder/${sample_id} \\
        ${params.fastqc_args} --memory $memory \\
        --threads ${task.cpus} \\
        ${read1} ${read2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """

}
