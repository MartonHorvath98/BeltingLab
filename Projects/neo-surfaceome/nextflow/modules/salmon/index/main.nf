process SALMON_INDEX {
    label 'thicc'
    tag "Salmon index on $genome"
    storeDir params.reference

    input:
    path transcriptome
    path genome

    output:
    path '*.salmon.idx', emit: index
    path 'versions.yml', emit: version

    script:
    def index_name = genome.name + ".salmon.idx"
    """
    ${module}
    
    # Create decoy file
    grep "^>" ${genome} | cut -d " " -f 1 | sed -e 's/>//g' > decoys.txt
    # Concatenate transcriptome and genome files
    cat ${transcriptome} ${genome} > gentrome.fa
    # Create Salmon index
    salmon index \\
        --transcripts gentrome.fa \\
        --decoys decoys.txt \\
        --index "${index_name}" \\
        -p $task.cpus \\
        ${params.salmon_index_args}

    rm -rf gentrome.fa decoys.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
