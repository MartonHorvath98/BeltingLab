process FASTQC {
    tag "FastQC on $sample_id"
    label "low_effort"

    publishDir "${params.report_outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(read1), path(read2)
    val subfolder

    output:
    path "$subfolder/${sample_id}/*_fastqc.zip"

    script:
    def memory = "${task.memory.toMega()}"
    """
    mkdir -p $subfolder/${sample_id}
    fastqc \\
        -o $subfolder/${sample_id} \\
        ${params.fastqc_args} --memory $memory \\
        --threads ${task.cpus} \\
        ${read1} ${read2}
    """

}
