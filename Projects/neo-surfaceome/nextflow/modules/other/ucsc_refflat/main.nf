process GTF2REFFLAT {
    tag "Convert GTF to genePred"
    label "low_effort"
    publishDir "${params.reference}", mode: 'copy', overwrite: true

    output:
    path "ref_annot.refFlat.txt.gz", emit: ref_flat
    path "versions.yml"		   , emit: version

    script:
    def VERSION="421" //Warning: UCSC does not provide version information on CLI. Update this string if needed!"
    """
    # Convert GTF to genePred format
    gtfToGenePred \\
        -genePredExt \\
        -geneNameAsName2 \\
        "${params.gtf_annot}" \\
        "ref_annot.refFlat.tmp.txt"
    
    # Parse genePred file to UCSC refFlat format
    paste <(cut -f 12 ref_annot.refFlat.tmp.txt) <(cut -f 1-10 ref_annot.refFlat.tmp.txt) > ref_annot.refFlat.txt
    rm ref_annot.refFlat.tmp.txt
    gzip ref_annot.refFlat.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}