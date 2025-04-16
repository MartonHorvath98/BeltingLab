process TRIM_GALORE {
    tag "Trim-galore on $sample"
    label "TRIM"

    publishDir params.trim_outdir, mode: 'copy'

    input:
    tuple val(sample), file(read1), file(read2)

    output:
    path '*_val_1.fq.gz', emit: trimmed_1
    path '*_val_2.fq.gz', emit: trimmed_2
    path '*_trimming_report.txt', emit: fastqc_report

    script:
    """
    echo "READ1: ${read1}"
    echo "READ2: ${read2}"
    
    trim_galore ${read1} ${read2} \\
     -j $task.cpus \\
     $params.trim_args \\
     --basename $sample
    """

}