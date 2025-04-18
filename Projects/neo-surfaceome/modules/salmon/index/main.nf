process INDEX {
    label "SALMON"
    tag "Salmon index on $genome"
    storeDir params.reference

    input:
    path transcriptome
    path genome

    output:
    path '${genome.simpleName}.salmon.idx', emit: index

    script:
    def index_name = genome.name + ".salmon.idx"
    """
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
        ${params.index_args}

    rm -rf gentrome.fa decoys.txt
    """
}
