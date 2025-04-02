process SALMON_INDEX {
    tag "salmon_index"    
    label "process_medium"

    conda "${moduleDir}/environment.yml"


    input:
    path(genome_fasta), path(transcript_fasta)

    output:
    path "${genome_fasta}salmon.idx/", emit: index

    when:
    file(genome_fasta).exists() && file(transcript_fasta).exists() && !file("${genome_fasta}salmon.idx/").exists()

    script:
    def args = task.ext.args ?: ''
    def decoys = ''
    def fasta = transcript_fasta
    if (genome_fasta){
        if (genome_fasta.endsWith('.gz')) {
            genome_fasta = "<(gunzip -c $genome_fasta)"
        }
        decoys='-d decoys.txt'
        gentrome='gentrome.fa'
    }
    if (transcript_fasta.endsWith('.gz')) {
        transcript_fasta = "<(gunzip -c $transcript_fasta)"
    }
    """
    if [ -n '$genome_fasta' ]; then
        grep '^>' $genome_fasta | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 | sed 's/>//g' > decoys.txt
        cat $transcript_fasta $genome_fasta > $gentrome
    fi

    salmon \\
        index \\
        --threads $task.cpus \\
        -t $gentrome \\
        $decoys \\
        $args \\
        -i $index
}