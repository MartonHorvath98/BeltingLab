process SEQKIT {
    tag "Seqkit on $sample_id"
    label "low_effort"

    publishDir "${params.report_outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(read1), path(read2)
    val subfolder

    output:
    path "$subfolder/${sample_id}_stats.tsv"

    script:
    """
    mkdir $subfolder
    seqkit stats ${read1} -j ${task.cpus} -Ta >> ${subfolder}/R1.txt
    seqkit stats ${read2} -j ${task.cpus} -Ta >> ${subfolder}/R2.txt

    cat  ${subfolder}/R1.txt  ${subfolder}/R2.txt |\\
        awk 'NR==1 || NR%2==0' |\\
        awk -F'\\t' -v OFS='\\t' '{sub(".*/", "", \$1)} 1' > $subfolder/${sample_id}_stats.tsv

    rm ${subfolder}/R1.txt
    rm ${subfolder}/R2.txt
    """

}
