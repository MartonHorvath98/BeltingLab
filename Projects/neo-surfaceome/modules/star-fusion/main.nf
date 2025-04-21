  process FUSION {
    tag "STAR-Fusion on $sample_id"
    label "STAR"

    publishDir params.fusion_outdir, mode: 'copy'

    input:
    tuple path(genome_lib_dir), val(sample_id), file(trim1), file(trim2), file(junction)

 
    output:
    path("star-fusion.fusion_predictions.*"), emit: fusion_preds

    script:
    """
    STAR-Fusion \\
        -J "${junction}" \\
        --left_fq "${trim1}" \\
        --right_fq "${trim2}" \\
        --genome_lib_dir "${genome_lib_dir}" \\
        --CPU $task.cores \\
        $params.fusion_args \\
        --output_dir "." --outTmpDir "tmp/"
    """
}