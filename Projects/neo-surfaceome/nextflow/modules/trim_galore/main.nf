process TRIM_GALORE {
    label 'low_effort'
    tag "Trim-galore on $sample"
    publishDir params.trim_outdir, mode: 'symlink'

    input:
    tuple val(sample), file(read1), file(read2)

    output:
    tuple val(sample), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
    tuple val(sample), path("*.txt"), emit: log
    tuple val(sample), path("*.zip"), emit: zip, optional: true
    tuple val(sample), path("*.html"), emit: html, optional: true
    path "*versions.yml", emit: versions


    script:
    """    
    trim_galore ${read1} ${read2} \\
     -j $task.cpus \\
     --gzip \\
     $params.trim_args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """

}