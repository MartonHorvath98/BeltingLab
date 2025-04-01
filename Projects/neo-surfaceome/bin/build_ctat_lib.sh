#!/bin/bash
# Using getopt

#####################################################
## Marton Horvath, JUN 2024 - last modified: MAR 2025
#####################################################
# ------------ Custom functions ------------ #
# Download function
get() {
	file=$1
	if ! wget --version >/dev/null 2>/dev/null ; then
		if ! curl --version >/dev/null 2>/dev/null ; then
			echo "Please install wget or curl somewhere in your PATH"
			exit 1
		fi
		curl -o $2 $1
		return $?
	else
		wget $1 -P $2
		return $?
	fi
}

# ------------- Load environment --------------- #
conda init bash
# Debugging: Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: Conda is not installed or not in PATH."
    exit 1
fi

# Debugging: Check if buildenv.yml exists
if [ ! -f "./config/buildenv.yml" ]; then
    echo "Error: ./config/buildenv.yml not found."
    exit 1
fi

# Attempt to create the environment
if [[ $(conda env list | grep "ctat-build") ]]; then
   conda activate ctat-build
else 
   conda env create --name ctat-build --file ./config/buildenv.yml
   conda activate ctat-build
   exit
fi

# ------------ Set global variables ------------ #
samtools=`which samtools`
prep_genome_lib=`which prep_genome_lib.pl`
star=`which STAR`
# Set the number of cores
threads=16
eval echo "Number of cores: ${threads}"
# Set the output directory
output="./ctat"

# DOWNLOAD REFERENCE GENOME AND ANNOTATION
echo "#############################################"
echo "# 1. Download reference genome & annotation #"
echo "#############################################"
# Set up the Ensembl release version
release=113
# Set up the Ensembl base URLs
ENSEMBL_GENOME=ftp://ftp.ensembl.org/pub/release-${release}/fasta/homo_sapiens/dna
ENSEMBL_TRANSCRIPTOME=ftp://ftp.ensembl.org/pub/release-${release}/fasta/homo_sapiens/cds
# Set up GENCODE annotation URLs
GENCODE_GFF=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/
GENCODE_GTF=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/
# Download the GRCh38 reference genome fasta file (unless it already exists)
GENOME=Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
if [ ! -f "${output}/Homo_sapiens.GRCh38.fa" ] ; then
	get ${ENSEMBL_GENOME}/$GENOME ${output} || (echo "Error getting $GENOME" && exit 1)
	gunzip ${output}/$GENOME || (echo "Error unzipping $GENOME" && exit 1)
	mv ${output}/${GENOME%.gz} ${output}/Homo_sapiens.GRCh38.fa
fi
GENOME="${output}/Homo_sapiens.GRCh38.fa"

# Download the GRCh38 reference genome fasta file (unless it already exists)
TRANSCRIPTOME=Homo_sapiens.GRCh38.cds.all.fa.gz
if [ ! -f "${output}/Homo_sapiens.GRCh38.cds.fa" ] ; then
	get ${ENSEMBL_TRANSCRIPTOME}/$TRANSCRIPTOME ${output} || (echo "Error getting $TRANSCRIPTOME" && exit 1)
	gunzip ${output}/$TRANSCRIPTOME || (echo "Error unzipping $TRANSCRIPTOME" && exit 1)
	mv ${output}/${TRANSCRIPTOME%.gz} ${output}/Homo_sapiens.GRCh38.cds.fa
fi
TRANSCRIPTOME="${output}/Homo_sapiens.GRCh38.cds.fa"

# Download the reference's gene annotation in GFF3 format (unless it already exists)
GFF="gencode.v47.primary_assembly.annotation.gff3.gz"
if [ ! -f "${output}/gencode.v47.gff" ] ; then
	get ${GENCODE_GFF}${GFF} ${output} || (echo "Error getting $GFF" && exit 1)
	gunzip ${output}/${GFF} || (echo "Error unzipping $GFF" && exit 1)
	mv ${output}/${GFF%.gz} ${output}/gencode.v47.gff
fi
GFF="${output}/gencode.v47.gff"

# Download the reference's gene annotation in GTF format (unless it already exists)
GTF="gencode.v47.primary_assembly.annotation.gtf.gz"
if [ ! -f "${output}/gencode.v47.gtf" ] ; then
	get ${GENCODE_GFF}${GTF} ${output} || (echo "Error getting $GTF" && exit 1)
	gunzip ${output}/${GTF} || (echo "Error unzipping $GTF" && exit 1)
	mv ${output}/${GTF%.gz} ${output}/gencode.v47.gtf
fi
GTF="${output}/gencode.v47.gtf"

# BUILD REFERENCE GENOME INDECES
echo "#############################################"
echo "# 2. Index reference genome with samtools   #"
echo "#############################################"
if "$samtools faidx ${GENOME}" ; then
	echo "Genome index built"
else
	echo "Index building failed; see error message"
fi

echo "###############################################"
echo "# 3. reate reference genom index              #"
echo "###############################################"
# Create STAR index directory
BUILD_INDEX="${prep_genome_lib} --genome_fa ${GENOME} --gtf ${GTF} --dfam_db human --pfam_db current --human_gencode_filter --CPU ${threads} --output_dir ${output}"
if "${BUILD_INDEX}" ; then
	echo "Genome index built"
else
	echo "Index building failed; see error message"
fi

echo "###############################################"
echo "# 4. Build STAR index                         #"
echo "###############################################"
# Create STAR index directory
STAR_DIR=${output}/ref_genome.star.idx/ && mkdir -p ${STAR_DIR}
# Build STAR index
STAR="${star} --runThreadN ${threads} --runMode genomeGenerate --genomeDir ${STAR_DIR}\
--genomeFastaFiles ${GENOME} --sjdbGTFfile ${GFF%.gff}.gtf --sjdbOverhang 100"
echo Running $STAR
if $STAR ; then
	echo "STAR index built"
else
	echo "Index building failed; see error message"
fi

       
# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit

