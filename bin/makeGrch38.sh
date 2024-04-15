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
## 4. Salmon
SALMON_EXE=./salmon
if [ ! -x "$SALMON_EXE" ] ; then
	if ! which salmon ; then
		echo "Could not find salmon in current directory or in PATH"
		exit 1
	else
		SALMON_EXE=`which salmon`
	fi
fi
## 5. Salmon
BED_EXE=./gff2bed
if [ ! -x "$BED_EXE" ] ; then
	if ! which gff2bed ; then
		echo "Could not find bedops in current directory or in PATH"
		exit 1
	else
		SALMON_EXE=`which gff2bed`
	fi
fi

# Download the GRCh38 reference genome fasta file (unless it already exists)
GENOME=Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
if [ ! -f "${dir}/${GENOME%.dna.primary_assembly.fa.gz}.fa" ] ; then
	get ${ENSEMBL_GENOME}/$GENOME ${dir} || (echo "Error getting $GENOME" && exit 1)
	gunzip ${dir}/$GENOME || (echo "Error unzipping $GENOME" && exit 1)
	mv ${dir}/${GENOME%.gz} ${dir}/${GENOME%.dna.primary_assembly.fa.gz}.fa
fi
GENOME=${dir}/${GENOME%.dna.primary_assembly.fa.gz}.fa

# Download the GRCh38 reference genome fasta file (unless it already exists)
TRANSCRIPTOME=Homo_sapiens.GRCh38.cds.all.fa.gz
if [ ! -f "${dir}/${TRANSCRIPTOME%.cds.all.fa.gz}.cds.fa" ] ; then
	get ${ENSEMBL_TRANSCRIPTOME}/$TRANSCRIPTOME ${dir} || (echo "Error getting $TRANSCRIPTOME" && exit 1)
	gunzip ${dir}/$TRANSCRIPTOME || (echo "Error unzipping $TRANSCRIPTOME" && exit 1)
	mv ${dir}/${TRANSCRIPTOME%.gz} ${dir}/${TRANSCRIPTOME%.cds.all.fa.gz}.cds.fa
fi
TRANSCRIPTOME=${dir}/${TRANSCRIPTOME%.cds.all.fa.gz}.cds.fa

# Download the reference's gene annotation in GFF3 format (unless it already exists)
GFF=Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gff3.gz
if [ ! -f "${dir}/${GFF%.${ENSEMBL_RELEASE}.gff3.gz}.gff" ] ; then
	get ${ENSEMBL_GFF3_BASE}/$GFF ${dir} || (echo "Error getting $GFF" && exit 1)
	gunzip ${dir}/$GFF || (echo "Error unzipping $GFF" && exit 1)
	mv ${dir}/${GFF%.gz} ${dir}/${GFF%.${ENSEMBL_RELEASE}.gff3.gz}.gff
fi
GFF=${dir}/${GFF%.${ENSEMBL_RELEASE}.gff3.gz}.gff

# BUILD REFERENCE GENOME INDECES
echo "##########################################################"
echo "# 1. Index reference genome with samtools                #"
echo "##########################################################"

SAMTOOLS="${SAMTOOLS_EXE} faidx ${GENOME}"
echo Running $SAMTOOLS
if $SAMTOOLS ; then
	echo "Genome index built"
else
	echo "Index building failed; see error message"
fi

echo "##########################################################"
echo "# 2. Convert GFF file to GTF and BED                     #"
echo "##########################################################"

GFFREAD="${GFFREAD_EXE} ${GFF} -T -o ${GFF%.gff}.gtf"
echo Running $GFFREAD
if $GFFREAD ; then
	gff2bed < ${GFF%.gff}.gtf > ${GFF%.gff}.bed
	echo "GFF to GTF conversion"
else
	echo "Conversion failed; see error message"
fi

echo "##########################################################"
echo "# 3. Extract splice sites and exons                      #"
echo "##########################################################"
# extract splice sites
HISAT_DIR=${dir}/hisat && mkdir -p ${HISAT_DIR}
hisat2_extract_splice_sites.py ${GFF%.gff}.gtf > ${HISAT_DIR}/${base_name}.ss
SS=${HISAT_DIR}/${base_name}.ss
# extract exons
hisat2_extract_exons.py ${GFF%.gff}.gtf > ${dir}/hisat/${base_name}.exon
EXON=${HISAT_DIR}/${base_name}.exon

echo "##########################################################"
echo " 4. Building HISAT2 index                                #"
echo "##########################################################"

HISAT="${HISAT2_BUILD_EXE} -p ${threads} --seed 100 ${GENOME} --ss ${SS} --exon ${EXON} ${HISAT_DIR}/${base_name}"
echo Running $HISAT
if $HISAT ; then
	echo "genome index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi

echo "##########################################################"
echo "# 5. Build Salmon index                                  #"
echo "##########################################################"
# Create decoys.txt file
SALMON_DIR=${dir}/salmon_index && mkdir -p ${SALMON_DIR}
grep "^>" ${GENOME} | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
# Concatenate transcriptome and genome
cat ${TRANSCRIPTOME} ${GENOME} > gentrome.fa
# Build Salmon index
SALMON_INDEX="${SALMON_EXE} index -t gentrome.fa -d decoys.txt -i ${SALMON_DIR} -k 13 -p ${threads}"
echo Running $SALMON_INDEX
if $SALMON_INDEX ; then
	rm gentrome.fa
	echo "Salmon index built"
else
	echo "Index building failed; see error message"
fi

# Remove concatenated transcriptome and genome
rm gentrome.fa

echo "##########################################################"
echo "# 6. Indexing complete                                   #"
echo "##########################################################"
echo "You may now use the reference genome with HISAT2, Salmon, and other tools."
