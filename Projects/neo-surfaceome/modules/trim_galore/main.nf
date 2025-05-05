process TRIM_GALORE {
    label 'low_effort'
    tag "Trim-galore on $sample"
    publishDir params.trim_outdir, mode: 'copy'

    input:
    tuple val(sample), file(read1), file(read2)

    output:
    tuple val(sample), file('*_val_1.fq.gz'), file('*_val_2.fq.gz'), emit: trimmed_reads
    path '*_trimming_report.txt', emit: fastqc_report

    script:
    def module = params.environment == 'uppmax' ? 'module load TrimGalore/0.6.1' : ''
    """
    ${module}
    echo "TRIM_GALORE run on executor ${task.executor}"
    echo "READ1: ${read1}"
    echo "READ2: ${read2}"
    name=\$(sed 's/Sample_//g' <<< ${sample})
    
    trim_galore ${read1} ${read2} \\
     -j $task.cpus \\
     $params.trim_args \\
     --basename \$name
    """

}