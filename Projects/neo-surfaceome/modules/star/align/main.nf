process ALIGNMENT {
    tag "STAR alignment on $sample_id"
    label "STAR"

    publishDir params.star_outdir, mode: 'copy'

    input:
    tuple path(star_index), val(sample_id), file(trim1), file(trim2)

 
    output:
    tuple val(sample_id), path("${sample_id}-Aligned.sortedByCoord.out.bam"), path("${sample_id}-Aligned.sortedByCoord.out.bai"), emit: out_bam
    path("${sample_id}-Chimeric.out.junction"), emit: junction
    path("${sample_id}-Log.final.out"), emit: log

    script:
    """
    STAR \\
	    --runThreadN "$task.cpus" \\
	    --genomeDir "${star_index}" \\
	    --readFilesIn ${trim1} ${trim2} \\
	    --outFileNamePrefix "${sample_id}-" \\
        ${params.star_args}

    samtools index -@ "$task.cpus" "${sample_id}-Aligned.sortedByCoord.out.bam"
    """
}
