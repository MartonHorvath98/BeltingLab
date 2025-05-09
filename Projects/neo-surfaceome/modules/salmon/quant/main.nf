process QUANTIFICATION {    
    label 'thicc'
    tag "Salmon on $sample_id"
    publishDir params.quant_outdir, mode: 'copy'

    input:
    tuple path(salmon_index), val(sample_id), file(trim1), file(trim2)

 
    output:
    tuple val(sample_id), path("${sample_id}"), emit: quant_files
    tuple val(sample_id), path("*info.json"), emit: json_info, optional: true
    tuple val(sample_id), path("*lib_format_counts.json"), emit: lib_format_counts, optional: true
    path "versions.yml", emit: versions
    
    
    script:
    """
    salmon quant --threads $task.cpus \\
        $params.quant_args \\
        -i $salmon_index \\
        -1 ${trim1} -2 ${trim2} \\
        -o $sample_id

    if [ -f ${sample_id}/aux_info/meta_info.json ]; then
        cp ${sample_id}/aux_info/meta_info.json "${sample_id}_meta_info.json"
    fi
    if [ -f ${sample_id}/lib_format_counts.json ]; then
        cp ${sample_id}/lib_format_counts.json "${sample_id}_lib_format_counts.json"
    fi

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}
