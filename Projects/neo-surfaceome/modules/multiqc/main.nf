process MULTIQC {
    label "MULTIQC"
    publishDir params.report_outdir, mode:'copy'

    input:
    path  multiqc_inputs

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}
