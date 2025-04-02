process HMM_SCAN {
    tag "$meta.id"
    label "process_medium"
    
    conda "${moduleDir}/environment.yml"

    input:
    path(hmm), path(query)

    output:
    path('*.tblout.dat'), emit: output
    path('*.log'), emit: log
    
    when:
    file(query).exists() && file(hmm).exists() && !file(output).exists()

    script:
    def args    = task.ext.args ?: args.default
    output_dir  = `dirname ${query}`
    query_name  = `basename ${query}`
    output   = "--domtblout ${output_dir}/${query_name}.tblout.dat"
    log    = "${output_dir}/${query_name}.log"
    """
    hmmscan \\
        $args \\
        --cpu $task.cpus \\
        $tbl_file \\
        $hmm \\
        $query \\
        2> $log \\
        | tee $log
    """
}