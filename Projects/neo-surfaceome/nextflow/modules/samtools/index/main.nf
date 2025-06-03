process SAMTOOLS_INDEX {
    tag "Samtools index on $sample_id"
    label "low_effort"

    publishDir "${params.star_outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(bam_input)

    output:
    tuple val(sample_id), path("*.bai") , optional:true, emit: bai
    path  "versions.yml"           , emit: versions

    script:
    """
    samtools \\
        index \\
        -@ ${task.cpus} \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}
