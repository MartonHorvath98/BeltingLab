process VEP_ANNOTATE {
    label 'thicc'
    tag "Annotate variants on $sample (VEP)"
    publishDir params.snv_outdir, mode: 'copy', overwrite: true

    input:
    path ref_genome
    tuple val(sample_id), file(vcf)

    output:
    tuple val(sample_id), path("*-vep.tsv")       , optional:true, emit: tab
    tuple val(sample_id), path("*-vep.html")      , optional:true, emit: report
    path "versions.yml"                           , emit: version


    script:
    """
    vep \\
        --input_file ${vcf} \\
        --format "vcf" \\
        --output_file ${sample_id}-vep.tsv \\
        --tab \\
	--stats_html --stats_file ${sample_id}-vep.html \\
	--offline --cache --dir_cache \$VEP_CACHE \\
	--fasta ${ref_genome} \\
	--assembly "GRCh38" \\
	--force \\
	--everything

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensembl-vep: \$(echo \$(vep --help 2>&1) | grep -oP "(?<=ensembl-vep : ).+?(?=\\ )")
    END_VERSIONS
    """

}

