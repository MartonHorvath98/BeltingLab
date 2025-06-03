process BCFTOOLS_STATS {
    label 'low_effort'
    tag "VC-stats on $sample (BCFtools)"
    publishDir params.snv_outdir, mode: 'copy', overwrite: true

    input:
    path ref_genome
    tuple val(sample), path(vcf)

    output:
    tuple val(sample), path("*stats.txt")  , emit: stats
    path  "versions.yml"                   , emit: version

    script:
    """
    bcftools stats \\
        --threads ${task.cpus} \\
        -F ${ref_genome} -s - \\
        ${vcf} > "${sample}.bcftools_stats.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}