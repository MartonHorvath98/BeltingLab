process STAR_INDEX {
    label "thicc"
    tag "STAR index on $genome"
    storeDir params.reference

    input:
    path genome
    path gtf

    output:
    path '${genome.simpleName}.star.idx', emit: index
    path 'versions.yml'			, emit: version

    script:
    def index_name = genome.name + ".star.idx"
    """
    STAR --runMode genomeGenerate \\
         --runThreadN $task.cpus \\
         --genomeDir $index_name \\
         --genomeFastaFiles $genome \\
         --sjdbGTFfile $gtf \\
         ${params.star_index_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR: \$(STAR --version | sed 's/STAR_//')
    END_VERSIONS
    """
}
