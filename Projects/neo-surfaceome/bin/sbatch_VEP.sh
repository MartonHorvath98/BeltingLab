#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p node -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J annotate-vcf
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4

# Marton Horvath, FEB 2025 - last modified: MAR 2025
# ------------ Start a screen environment ----- #
screen -S VEPAnnotation
# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools vep/113.0

# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/

# ------------ Set global variables ------------ #
# Path to the variant calling result directory
SAMPLES="$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')"
# Fasta reference
REF="${WHARF}refs/ctat/"
# Input directory
INPUT_DIR="${WHARF}results/05_SNV/"
eval echo "Variant files from: ${INPUT_DIR}"
# Vep cache directory
CACHE="${VEP_CACHE}"
eval echo "VEP CACHE is in: ${CACHE}"
# Vep custom annotation
VEP_RESOURCE="${WHARF}resources/vep/"
eval echo "VEP custom annotation: ${VEP_RESOURCE}"

# ----------- Check executables ---------------- #
## 1. VEP 
VEP_EXE=./vep
if [ ! -x "$VEP_EXE" ] ; then
	if ! which vep; then
		echo "Could not find ENSEMBL's variant effect predictor (VEP) in current directory or in PATH"
		exit 1
	else
		VEP_EXE=`which vep`
	fi
fi

echo "#################################################"
echo "# Single nucleotide variant annotation with VEP #"
echo "#################################################"
for sample_dir in $(find ${INPUT_DIR} -maxdepth 1 -mindepth 1 -type d)
do
    variant_caller=$(basename ${sample_dir})
    
    OUTPUT_DIR="${INPUT_DIR}${variant_caller}/"
    if [ ! -d ${OUTPUT_DIR} ] ; then
	mkdir -p ${OUTPUT_DIR}
    fi

    for sample in ${SAMPLES}
    do
	echo -e "Annotation variants called by ${variant_caller} on: ${sample}"
	srun --exclusive --ntasks=4 --cpus-per-task=4 "${VEP_EXE}" --input_file "${sample_dir}/${sample}-${variant_caller}.vcf.gz" \
	    --format "vcf" --output_file "${OUTPUT_DIR}${sample}-${variant_caller}-VEP.tsv" --tab \
	    --stats_html --stats_file "${OUTPUT_DIR}${sample}-${variant_caller}-VEP-summary.html" \
	    --offline --cache --dir_cache "${CACHE}" \
	    --fasta "${REF}ref_genome.fa" \
	    --assembly "GRCh38" \
	    --fork 4 \
	    --force \
	    --everything
	    
    done
done


	    # --gtf "${REF}ref_annot.gtf.gz" \
	    #--distance 150 \ 
	    # --check_existing \
	    # --canonical \
	    # --sift b \
	    # --polyphen b \
	    # --regulatory \
	    # --domains \
	    # --gene_phenotype \
	    # --symbol \
	    # --protein \
	    # --uniprot \
	    # --af --af_gnomade \

	    # --custom file=${VEP_RESOURCE}hg38.phyloP100way.bw,short_name=PhyloP,format=bigwig,type=exact \
	    # --custom file=${VEP_RESOURCE}hg38.phastCons100way.bw,short_name=PhastCons,format=bigwig,type=exact \
	    # --custom file=${VEP_RESOURCE}encode_rna_binding_sorted.bed.gz,short_name=Encode,format=bed,type=overlap \
	    # --custom file=${VEP_RESOURCE}eve_data.vcf.gz,short_name=EVE,format=vcf,fields="EVE%ProtMut%Class10%Class50%Class90" \
	    # --custom file=${VEP_RESOURCE}hg38.phyloP100way.bw,short_name=PhyloP,format=bigwig,type=exact \
	    # --custom file=${VEP_RESOURCE}uniprot_human_domain_reformatted.bed.gz,short_name=UniprotDom,format=bed,type=overlap \
	    # --custom file=${VEP_RESOURCE}uniprot_human_topo_dom_reformatted.bed.gz,short_name=UniprotTopo,format=bed,type=overlap \
	    # --custom file=${VEP_RESOURCE}clinvar_20240813.vcf.gz,short_name=ClinVar,format=vcf,fields="ALLELEID%CLNDISDB% \
	    # CLNDN%CLNHGVS%CLNSIG%CLNREVSTAT%CLNSIGCONF%CLNSIGINCL%CLNVC%CLNVCSO%CLNVI%DBVARID%GENEINFO%MC%ORIGIN" \
	    # --plugin MaxEntScan,${VEP_RESOURCE}MaxEnTScan,SWA,NCSS \
	    # --plugin AlphaMissense,${VEP_RESOURCE}AlphaMissense_hg38.tsv.gz \
	    # --plugin SpliceVault,${VEP_RESOURCE}SpliceVault_GRCh38_data.tsv.gz \
# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
		
	    






