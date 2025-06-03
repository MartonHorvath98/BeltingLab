process STAR_ALIGN {
    tag "STAR alignment on $sample"
    label "thicc"

    publishDir params.star_outdir, mode: 'symlink'

    input:
    tuple path(star_index), val(sample), file(trim1), file(trim2)

 
    output:
    tuple val(sample), path("*-Aligned.sortedByCoord.out.bam"), emit: out_bam
    tuple val(sample), path("*-Chimeric.out.junction")        , emit: junction
    tuple val(sample), path("*-Log.final.out")                , emit: log
    path "versions.yml"                                          , emit: version

    script:
    """
    STAR \\
		--runThreadN "$task.cpus" \\
		--genomeDir "${star_index}" \\
		--readFilesCommand "gunzip -c" \\
		--readFilesIn ${trim1} ${trim2} \\
		--outFileNamePrefix "${sample}-" \\
		${star_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR: \$(STAR --version | sed 's/STAR_//')
    END_VERSIONS
    """
}
