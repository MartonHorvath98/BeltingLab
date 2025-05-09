process AGAT_CONVERT {
    publishDir params.reference, mode: 'copy'

    output:
    path "ref_annot.gff", emit: gff_file
    path "ref_annot.bed", emit: bed_file

    script:
    """
    #STEP 1: convert gtf annotation to gff3 format
    agat_convert_sp_gxf2gxf.pl -g ${params.gtf_annot} -o ref_annot.gff
    
    #STEP 2: convert gff3 annotation to bed format
    agat_convert_sp_gxf2bed.pl -g ref_annot.gff -o ref_annot.bed
    """
}