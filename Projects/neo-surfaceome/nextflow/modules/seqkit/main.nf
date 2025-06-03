process SEQKIT {
    tag "Seqkit on $sample"
    label "low_effort"

    publishDir "${params.report_outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)
    val subfolder

    output:
    tuple val(sample), path("$subfolder/${sample}_stats.tsv"), emit: stats
    path "versions.yml"	  		                     , emit: version

    script:
    """
    mkdir $subfolder
    seqkit stats ${read1} -j ${task.cpus} -Ta >> ${subfolder}/R1.txt
    seqkit stats ${read2} -j ${task.cpus} -Ta >> ${subfolder}/R2.txt

    cat  ${subfolder}/R1.txt  ${subfolder}/R2.txt |\\
        awk 'NR==1 || NR%2==0' |\\
        awk -F'\\t' -v OFS='\\t' '{sub(".*/", "", \$1)} 1' > ${subfolder}/${sample}_stats.tsv

    rm ${subfolder}/R1.txt
    rm ${subfolder}/R2.txt

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """

}
