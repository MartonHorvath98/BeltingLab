process MULTIQC {
    label 'low_effort'
    publishDir params.report_outdir, mode:'copy'

    input:
    path  multiqc_inputs

    output:
    path 'multiqc_report.html'

    script:
    def module = params.environment == 'uppmax' ? 'module load MultiQC/1.10.1' : ''
    """
    ${module}
    multiqc .
    """
}
