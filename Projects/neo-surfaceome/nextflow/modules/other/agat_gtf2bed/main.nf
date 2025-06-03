process AGAT_CONVERT {
    input:
    path gtf_annot

    output:
    path "${ref_path}/Homo_sapiens.GRCh38.bed", emit: genes_bed

    script:
    def ref_path = gtf_annot.

    """
    #STEP 1: convert gtf annotation to gff3 format
    agat_convert_sp_gxf2gxf.pl -g ${ref_path}/Homo_sapiens.GRCh38.gff --out ${ref_path}/Homo_sapiens.GRCh38.bed
    fi
    """
}