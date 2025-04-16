process INDEX {
    label "SALMON"
    tag "Salmon index on $genome"
    publishDir params.reference, mode: 'copy'

    input:
    path transcriptome
    path genome

    output:
    path '${genome}.salmon.idx', emit: index

    script:
    """
    # mkdir -p ${genome}.salmon.idx
    # Create decoy file
    grep "^>" ${genome} | cut -d " " -f 1 | sed -e 's/>//g' > decoys.txt
    # Concatenate transcriptome and genome files
    cat ${transcriptome} ${genome} > gentrome.fa
    # Create Salmon index
    salmon index \\
        --transcripts gentrome.fa \\
        --decoys decoys.txt \\
        --index "${genome}.salmon.idx" \\
        -p $task.cpus \\
        ${params.index_args}

    rm -rf gentrome.fa decoys.txt
    """
}
