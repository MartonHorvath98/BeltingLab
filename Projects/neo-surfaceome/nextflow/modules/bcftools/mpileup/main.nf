process BCFTOOLS_MPILEUP {
    label 'thicc'
    tag "Variant calling on $sample (BCFtools)"
    publishDir params.snv_outdir, mode: 'copy', overwrite: true

    input:
    path ref_genome
    tuple val(sample), path(bam), path(bam_sub1), path(bam_sub2), path(bam_sub3)

    output:
    tuple val(sample), path("*vcf.gz")     , emit: vcf
    tuple val(sample), path("*vcf.gz.tbi") , emit: vcf_index
    path  "versions.yml"                   , emit: version

    script:
    """
    bcftools mpileup \\
        --threads "${task.cpus}" \\
        -Ou --fasta-ref "${ref_genome}" \\
        ${bam} ${bam_sub1} ${bam_sub2} ${bam_sub3} \\
        ${params.snv_args_mpileup} | bcftools call \\
        -m --variants-only -Oz > "${sample}.vcf.gz"

    tabix -p vcf -f ${sample}.vcf.gz
    rm ${bam} ${bam_sub1} ${bam_sub2} ${bam_sub3}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) |  grep -oP "(?<=Version: ).+?(?=\\ )")	
    END_VERSIONS
    """
}