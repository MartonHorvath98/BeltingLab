 process PICARD_DEDUP {
    label "generic_task"
    tag "Removing duplicates of ${sample}"
    publishDir "${params.snv_outdir}", mode: 'symlink'

    input:
    tuple val(sample), path(bam_file)

    output:
    tuple val(sample), path("*dedup.bam")       , emit: bam
    tuple val(sample), path("*-dup-metrics.txt"), emit: metrics
    path "versions.yml"                            , emit: version

    script:
    """
    # Collect RNA-seq metrics
    java -jar \$PICARD MarkDuplicates \\
		--INPUT ${bam_file} \\
		--OUTPUT ${sample}-dedup.bam \\
		--METRICS_FILE ${sample}-dup-metrics.txt \\
		--ASSUME_SORT_ORDER coordinate
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
	picard: \$(echo \$(java -jar \$PICARD MarkDuplicates 2>&1) | grep -oP "(?<=Version:).+?(?=\\ )")
    END_VERSIONS
    """
}