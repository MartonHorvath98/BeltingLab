process EXTRACT_RRNA {
    label "low_effort"
    tag "Extract rRNA intervals from reference genome"
    publishDir params.reference, mode: 'copy'


    output:
    path "rRNA_intervals.txt", emit: rrna_intervals
    path "versions.yml"	     , emit: version

    script:
    """
    # Sequence names and lengths. (Must be tab-delimited.)
    cut -f1,2 "${params.genome_file}.fai" | \\
        perl -lane 'print "\\@SQ\\tSN:\$F[0]\\tLN:\$F[1]\\tAS:GRCh38"' | \\
        grep -v _ \\
        >> "rRNA_intervals.txt"
    
    # Intervals for rRNA transcripts.
    grep 'gene_type "rRNA"' "${params.gtf_annot}" | \\
        awk '\$3 == "transcript"' | \\
        cut -f1,4,5,7,9 | \\
        perl -lane '
            /transcript_id "([^"]+)"/ or die "no transcript_id on \$.";
            print join "\\t", (@F[0,1,2,3], \$1)
        ' | \\
        sort -k1V -k2n -k3n \\
        >> "rRNA_intervals.txt"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | cut -d' ' -f3)
	perl: \$(perl --version | grep "This is perl" | grep -oP "(?<=\\(v).+?(?=\\))")
    END_VERSIONS
    """
}