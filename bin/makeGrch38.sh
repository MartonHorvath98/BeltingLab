#!/bin/bash

# Downloads sequence for the GRCh38 release version v.108 (2022 Oct)
# of H. sapiens (human) from Ensembl.
#
# Note that Ensembl's GRCh38 build has three categories of compressed fasta
# files:
# 'dna' - unmasked genomic DNA sequences.
# 'dna_rm' - masked genomic DNA.  Interspersed repeats and low complexity regions
#			 are detected with the RepeatMasker tool and masked with 'N's.
# 'dna_sm' - soft-masked genomic DNA. All repeats and low complexity regions have
#			 been replaced with lowercased versions of their nucleic base.
#
# By default, this script builds and index for just the unmasked files at the primary 
# assembly level, that contains all toplevel sequence regions excluding haplotypes
# and patches.
#! bin/bash

# Take user arguments: input directory
while getopts d:j:b: flag
do
    case "${flag}" in
        d) dir=${OPTARG};;
		j) threads=${OPTARG};;
		b) base_name=${OPTARG};;
    esac
done

# Set up the Ensembl release version and base URLs
ENSEMBL_RELEASE=108
ENSEMBL_FASTA_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna
ENSEMBL_GFF3_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gff3/homo_sapiens/

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

# Check command executable paths
## 1. Hisat2
HISAT2_BUILD_EXE=./hisat2-build
if [ ! -x "$HISAT2_BUILD_EXE" ] ; then
	if ! which hisat2-build ; then
		echo "Could not find hisat2-build in current directory or in PATH"
		exit 1
	else
		HISAT2_BUILD_EXE=`which hisat2-build`
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

# Download the GRCh38 reference genome fasta file (unless it already exists)
FASTA=Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
if [ ! -f "${dir}/${FASTA%.dna.primary_assembly.fa.gz}.fa" ] ; then
	get ${ENSEMBL_FASTA_BASE}/$FASTA ${dir} || (echo "Error getting $FASTA" && exit 1)
	gunzip ${dir}/$FASTA || (echo "Error unzipping $FASTA" && exit 1)
	mv ${dir}/${FASTA%.gz} ${dir}/${FASTA%.dna.primary_assembly.fa.gz}.fa
fi
FASTA=${dir}/${FASTA%.dna.primary_assembly.fa.gz}.fa

# Download the reference's gene annotation in GFF3 format (unless it already exists)
GFF=Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gff3.gz
if [ ! -f "${dir}/${GFF%.gff.gz}.gtf" ] ; then
	get ${ENSEMBL_GFF3_BASE}/$GFF ${dir} || (echo "Error getting $GFF" && exit 1)
	gunzip ${dir}/$GFF || (echo "Error unzipping $GFF" && exit 1)
	mv ${dir}/${GFF%.gz} ${dir}/${GFF%.${ENSEMBL_RELEASE}.gff3.gz}.gff
fi
GFF=${dir}/${GFF%.${ENSEMBL_RELEASE}.gff3.gz}.gff

# BUILD REFERENCE GENOME INDECES
## 1. Index reference genome with samtools
SAMTOOLS="${SAMTOOLS_EXE} faidx ${FASTA}"
if [ ! -f "${FASTA}.fai"] ; then
	echo Running $SAMTOOLS
	if $SAMTOOLS ; then
		echo "Genome index built"
	else
		echo "Index building failed; see error message"
	fi
fi
## 2. Convert GFF file to GTF
GFFREAD="${GFFREAD_EXE} ${GFF} -T -o ${GFF%.gff}.gtf"

if [ ! -f "${GFF%.gff}.gtf"] ; then
	echo Running $GFFREAD
	if $GFFREAD ; then
		echo "GFF to GTF conversion"
	else
		echo "Conversion failed; see error message"
	fi
fi
## 3. Build HISAT2 index
HISAT="${HISAT2_BUILD_EXE} -p ${threads} --seed 100 ${FASTA} ${dir}/${base_name}"
if [ ! -f "${dir}/${base_name}.1.ht2" ] ; then
	echo Running $HISAT
	if $HISAT ; then
		echo "genome index built; you may remove fasta files"
	else
		echo "Index building failed; see error message"
	fi
else
	echo "Index already exists; you may remove fasta files"
fi
