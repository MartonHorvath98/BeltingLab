#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p node -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J align-STAR

# Marton Horvath, JUN 2024 - last modified: DEC 2024
# ------------ Start a screen environment ----- #
screen -S STARalignment
# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools star/2.7.8a star-fusion/1.10.1 hmmer/3.3.2 blast/2.15.0+

# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/

CTAT_DIR=${WHARF}refs/ctat/
eval echo "Indexes are built: ${CTAT_DIR}"
# ------------ Set global variables ------------ #
# Set the number of cores
CORES=16
eval echo "Number of cores: ${CORES}"

# Set up the Ensembl release version and base URLs
ENSEMBL_RELEASE=108
ENSEMBL_GENOME=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna
ENSEMBL_GFF3_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gff3/homo_sapiens/
ENSEMBL_TRANSCRIPTOME=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/cds

# Custom download function
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
# ------- Check command executable paths -------
## 1. Genome build
BUILD_EXE=/sw/bioinfo/star-fusion/1.10.1/bianca/ctat-genome-lib-builder/prep_genome_lib.pl
if [ ! -x "$BUILD_EXE" ] ; then
    if ! which "$BUILD_EXE" ; then
        echo "Could not find utility script 'prep_genome_lip.pl' in current directory or in PATH"
        exit 1
    else
        BUILD_EXE=`which "$BUILD_EXE"`
    fi
fi

## 2. Samtools
SAMTOOLS_EXE=./samtools
if [ ! -x "$SAMTOOLS_EXE" ] ; then
	if ! which samtools ; then
		echo "Could not find samtools in current directory or in PATH"
		exit 1
	else
		SAMTOOLS_EXE=`which samtools`
	fi
fi

## 3. GFFread
GFFREAD_EXE=./gffread
if [ ! -x "$GFFREAD_EXE" ] ; then
	if ! which gffread ; then
		echo "Could not find gffread in current directory or in PATH"
		exit 1
	else
		GFFREAD_EXE=`which gffread`
	fi
fi

## 4. STAR
STAR_EXE=./STAR
if [ ! -x "$STAR_EXE" ] ; then
	if ! which STAR ; then
		echo "Could not find STAR in current directory or in PATH"
		exit 1
	else
		STAR_EXE=`which STAR`
	fi
fi


# DOWNLOAD REFERENCE GENOME AND ANNOTATION
echo "#############################################"
echo "# 1. Download reference genome & annotation #"
echo "#############################################"

# Download the GRCh38 reference genome fasta file (unless it already exists)
GENOME=Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
if [ ! -f "${CTAT_DIR}/${GENOME%.dna.primary_assembly.fa.gz}.fa" ] ; then
	get ${ENSEMBL_GENOME}/$GENOME ${CTAT_DIR} || (echo "Error getting $GENOME" && exit 1)
	gunzip ${CTAT_DIR}/$GENOME || (echo "Error unzipping $GENOME" && exit 1)
	mv ${CTAT_DIR}/${GENOME%.gz} ${CTAT_DIR}/${GENOME%.dna.primary_assembly.fa.gz}.fa
fi
GENOME=${CTAT_DIR}/${GENOME%.dna.primary_assembly.fa.gz}.fa

# Download the GRCh38 reference genome fasta file (unless it already exists)
TRANSCRIPTOME=Homo_sapiens.GRCh38.cds.all.fa.gz
if [ ! -f "${CTAT_DIR}/${TRANSCRIPTOME%.cds.all.fa.gz}.cds.fa" ] ; then
	get ${ENSEMBL_TRANSCRIPTOME}/$TRANSCRIPTOME ${CTAT_DIR} || (echo "Error getting $TRANSCRIPTOME" && exit 1)
	gunzip ${CTAT_DIR}/$TRANSCRIPTOME || (echo "Error unzipping $TRANSCRIPTOME" && exit 1)
	mv ${CTAT_DIR}/${TRANSCRIPTOME%.gz} ${CTAT_DIR}/${TRANSCRIPTOME%.cds.all.fa.gz}.cds.fa
fi
TRANSCRIPTOME=${CTAT_DIR}/${TRANSCRIPTOME%.cds.all.fa.gz}.cds.fa

# Download the reference's gene annotation in GFF3 format (unless it already exists)
GFF=Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gff3.gz
if [ ! -f "${CTAT_DIR}/${GFF%.${ENSEMBL_RELEASE}.gff3.gz}.gff" ] ; then
	get ${ENSEMBL_GFF3_BASE}/$GFF ${CTAT_DIR} || (echo "Error getting $GFF" && exit 1)
	gunzip ${CTAT_DIR}/$GFF || (echo "Error unzipping $GFF" && exit 1)
	mv ${CTAT_DIR}/${GFF%.gz} ${CTAT_DIR}/${GFF%.${ENSEMBL_RELEASE}.gff3.gz}.gff
fi
GFF=${CTAT_DIR}/${GFF%.${ENSEMBL_RELEASE}.gff3.gz}.gff

# BUILD REFERENCE GENOME INDECES
echo "#############################################"
echo "# 2. Index reference genome with samtools   #"
echo "#############################################"

SAMTOOLS="${SAMTOOLS_EXE} faidx ${GENOME}"
echo Running $SAMTOOLS
if $SAMTOOLS ; then
	echo "Genome index built"
else
	echo "Index building failed; see error message"
fi

echo "##############################################"
echo "# 3. Convert GFF file to GTF and BED         #"
echo "##############################################"

GFFREAD="${GFFREAD_EXE} ${GFF} -T -o ${GFF%.gff}.gtf"
echo Running $GFFREAD
if $GFFREAD ; then
	gff2bed < ${GFF%.gff}.gtf > ${GFF%.gff}.bed
	echo "GFF to GTF conversion"
else
	echo "Conversion failed; see error message"
fi

echo "###############################################"
echo "# 4. Build STAR index                         #"
echo "###############################################"
# Create STAR index directory
STAR_DIR=${CTAT_DIR}/star_index && mkdir -p ${STAR_DIR}
# Build STAR index
STAR="${STAR_EXE} --runThreadN ${threads} --runMode genomeGenerate --genomeDir ${STAR_DIR}\
--genomeFastaFiles ${GENOME} --sjdbGTFfile ${GFF%.gff}.gtf --sjdbOverhang 100"
echo Running $STAR
if $STAR ; then
	echo "STAR index built"
else
	echo "Index building failed; see error message"
fi

# ------------ Create reference genom index ---- #
# Create STAR index directory
BUILD_INDEX="${BUILD_EXE} --genome_fa  ${GENOME}  --gtf ${GFF%.gff}.gtf --dfam_db human --fusion_annot_lib ${CTAT_DIR}fusion_lib.Mar2021.dat.gz --human_gencode_filter --pfam_db current --CPU 16 --output_dir ${CTAT_DIR}"
if [ -d "$REF" ] 
then
    echo "Genome index $(basename $CTAT_DIR) has already been created"
else
    echo Running $BUILD_INDEX
    if $BUILD_INDEX ; then
	echo "Genome index built"
    else
	echo "Index building failed; see error message"
    fi
fi
       
# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit

