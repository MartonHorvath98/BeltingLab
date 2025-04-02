process HMM_PRESS {
    tag "$meta.id"
    label "process_single"
    
    conda "${moduleDir}/environment.yml"

    input:
    path(hmm)

    output:
    path("${hmm}.h3*"), emit: output
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    hmmpress \\
        ${hmm}
    """
}