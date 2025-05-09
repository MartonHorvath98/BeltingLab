process MULTIQC {
    label 'low_effort'
    publishDir params.report_outdir, mode:'copy'

    input:
    path  multiqc_files, stageAs: "?/*"

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    script:
    """
    multiqc \\
        --force \\
        .
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
