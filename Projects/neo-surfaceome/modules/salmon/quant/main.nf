process QUANTIFICATION {    
    label 'thicc'
    tag "Salmon on $sample_id"
    publishDir params.quant_outdir, mode: 'copy'

    input:
    tuple path(salmon_index), val(sample_id), file(trim1), file(trim2)

 
    output:
    tuple val(sample_id), path("${sample_id}"), emit: quant_files
    
    script:
    def module = params.environment == 'uppmax' ? 'module load Salmon/1.10.1' : ''
    """
    ${module}
    
    salmon quant --threads $task.cpus \\
        $params.quant_args \\
        -i $salmon_index \\
        -1 ${trim1} -2 ${trim2} \\
        -o $sample_id
    """
}
