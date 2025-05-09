process EXTRACT_RRNA {
    publishDir params.reference, mode: 'copy'

    output:
    path "rRNA_intervals.bed", emit: rrna_intervals

    script:
    """
    # Sequence names and lengths. (Must be tab-delimited.)
    cut -f1,2 "${params.genome_file}.fai" | \\
        perl -lane 'print "\\@SQ\\tSN:\$F[0]\\tLN:\$F[1]\\tAS:GRCh38"' | \\
        grep -v _ \\
        >> "rRNA_intervals.bed"
    
    # Intervals for rRNA transcripts.
    grep 'gene_type "rRNA"' "${params.gtf_annot}" | \\
        awk '\$3 == "transcript"' | \\
        cut -f1,4,5,7,9 | \\
        perl -lane '
            /transcript_id "([^"]+)"/ or die "no transcript_id on \$.";
            print join "\\t", (@F[0,1,2,3], \$1)
        ' | \\
        sort -k1V -k2n -k3n \\
        >> "rRNA_intervals.bed"
     
    """
}