 process PICARD_DOWNSAMPLE {
    label "generic_task"
    tag "Subsample ${sample}"
    publishDir "${params.snv_outdir}", mode: 'symlink'

    input:
    tuple val(sample), path(bam_file)

    output:
    tuple val(sample), path("*25pc.bam"), path("*50pc.bam"), path("*75pc.bam"), emit: sub_bams
    path "versions.yml"                                                          , emit: version

    script:
    """
    # Create a 25% subsample of the input BAM file
    java -jar \$PICARD DownsampleSam \\
      --INPUT ${bam_file} \\
      --OUTPUT "${sample}-25pc.bam" \\
      --RANDOM_SEED 42 \\
      --PROBABILITY 0.25 \\
      --VALIDATION_STRINGENCY SILENT

    # Create a 50% subsample of the input BAM file
    java -jar \$PICARD DownsampleSam \\
      --INPUT ${bam_file} \\
      --OUTPUT "${sample}-50pc.bam" \\
      --RANDOM_SEED 42 \\
      --PROBABILITY 0.50 \\
      --VALIDATION_STRINGENCY SILENT

    # Create a 75% subsample of the input BAM file
    java -jar \$PICARD DownsampleSam \\
      --INPUT ${bam_file} \\
      --OUTPUT "${sample}-75pc.bam" \\
      --RANDOM_SEED 42 \\
      --PROBABILITY 0.75 \\
      --VALIDATION_STRINGENCY SILENT
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(java -jar \$PICARD DownsampleSam 2>&1) | grep -oP "(?<=Version:).+?(?=\\ )")
    END_VERSIONS
    """
}
