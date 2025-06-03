  process FUSION {
    tag "STAR-Fusion on $sample"
    label "generic_task"

    publishDir params.fusion_outdir, mode: 'copy'

    input:
    path genome_lib_dir
    tuple val(sample), file(trim1), file(trim2), file(junction)

    output:
    tuple val(sample), path("*.fusion_predictions.tsv")                   , emit: fusions
    tuple val(sample), path("*.abridged.tsv")                             , emit: abridged
    tuple val(sample), path("*.coding_effect.tsv")     , optional: true   , emit: coding_effect
    path "versions.yml"                                                      , emit: version

    script:
    """
    STAR-Fusion \\
        --genome_lib_dir ${genome_lib_dir} \\
        --left_fq ${trim1} \\
        --right_fq ${trim2} \\
        -J ${junction} \\
        --CPU ${task.cpus} \\
        --output_dir . \\
        ${params.fusion_args}
    
    mv star-fusion.fusion_predictions.tsv ${sample}.starfusion.fusion_predictions.tsv
    mv star-fusion.fusion_predictions.abridged.tsv ${sample}.starfusion.abridged.tsv
    mv star-fusion.fusion_predictions.abridged.coding_effect.tsv ${sample}.starfusion.abridged.coding_effect.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """
}