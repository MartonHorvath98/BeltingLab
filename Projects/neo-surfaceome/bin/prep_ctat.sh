#!/bin/bash
# Using getopt


######################################################
# Marton Horvath, JUN 2024 - last modified: APR 2025 #
######################################################
# ------------ Custom functions ------------ #
# Usage function
usage() {
    echo "Usage: $0 -o <output_dir> -p <pfam_hmm>"
    exit 1
}
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
# --------------- Main script --------------- #
# Parse command-line options
while getopts "o:p:" opt; do
    case "$opt" in
		o) output=$OPTARG ;;
        p) pfam=$OPTARG ;;
        *) usage ;;
    esac
done
# QC for arguments 
if [[ -z "$pfam" || ! -f "$pfam" ]]; then
    echo "Error: Invalid or missing path (-p). Path should point to the Pfam-A.hmm file!"
    usage
fi

# Create and activate the environment
mamba env create -f config/buildenv.yml

# ------- Check command executable paths -------
## 1. Salmon index
SALMON_EXE=salmon
if [ ! -x "$SALMON_EXE" ] ; then
    if ! which "$SALMON_EXE" ; then
        echo "Could not find salmon in current directory or in PATH"
        exit 1
    else
        SALMON_EXE=`which salmon`
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
## 3. HMMSCAN
HMM_EXE=./hmmscan
if [ ! -x "$HMM_EXE" ] ; then
	if ! which hmmscan; then
		echo "Could not find 'hmmscan' in current directory or in PATH"
		exit 1
	else
		HMM_EXE=`which hmmscan`
	fi
fi


# DOWNLOAD THE CTAT GENOME LIBRARY
echo "#################################################"
echo "# 1. Download CTAT genome library (plug-n-play) #"
echo "#################################################"
# Set up the URL to the GRCh38_gencode_v37_CTAT_lib_Mar012021 plug-n-play library
ctat_url="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/"
# Create directory structure for Ensembl protein FASTA files
output_dir="${output}"
mkdir -p ${output_dir}

# Download the CTAT genome library
CTAT_LIB="GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz"
if [ ! -d "${output_dir}ctat_genome_lib/" ] ; then
	get ${ctat_url}${CTAT_LIB} ${output_dir} || (echo "Error getting ${CTAT_LIB}" && exit 1)
	tar -xvzf ${output_dir}${CTAT_LIB} || (echo "Error unzipping ${CTAT_LIB}" && exit 1)
    rm ${output_dir}${CTAT_LIB} || (echo "Error removing ${CTAT_LIB}" && exit 1)
    mv ${output_dir}${CTAT_LIB%.tar.gz}/ctat_genome_lib_build_dir/ ${output_dir}ctat_genome_lib/
    rm -r ${output_dir}${CTAT_LIB%.tar.gz} || (echo "Error removing ${CTAT_LIB%.tar.gz}" && exit 1)
else
	echo "${CTAT_LIB%.tar.gz} have been dowloaded already..."
fi
CTAT_LIB="${output_dir}ctat_genome_lib/"
echo -e "Downloaded CTAT plug-n-play genome library: ${CTAT_LIB}"

SALMON_IDX="${CTAT_DIR}ref_genome.fa.salmon.idx/"

# GENERATE INDEX FOR SALMON
echo "##########################################################"
echo "# 2. Build Salmon index                                  #"
echo "##########################################################"

GENOME="${CTAT_LIB}ref_genome.fa"
TRANSCRIPTOME="${CTAT_LIB}ref_annot.cds"
GTF="${CTAT_LIB}ref_annot.gtf"
# Create Salmon index directory
BUILD_INDEX="${SALMON_EXE} index -t ${CTAT_LIB}gentrome.fa -d ${CTAT_LIB}decoys.txt -i ${SALMON_IDX} -k 31 --gencode"

if [ ! -d "${RSALMON_INDEX}" ] ; then
    # Create decoys.txt file
    mkdir -p "${SALMON_IDX}"
    grep "^>" ${GENOME} | cut -d " " -f 1 > "${CTAT_LIB}decoys.txt"
    sed -i.bak -e 's/>//g' "${CTAT_LIB}decoys.txt"
    
    # Concatenate transcriptome and genome
    cat ${TRANSCRIPTOME} ${GENOME} > "${CTAT_LIB}gentrome.fa"

    # Build Salmon index
    echo "Running $BUILD_INDEX"
    if $BUILD_INDEX ; then
	    # Remove concatenated transcriptome and genome
	    rm -f "${CTAT_LIB}decoys.txt.bak" "${CTAT_LIB}gentrome.fa"
	    echo "Salmon index built"
    else
	    echo "Index building failed; see error message"
    fi
else
    echo "Salmon index has already been built..."
fi
       
# ANNOTATE WITH PFAM USING HMMSCAN
echo "############################################################"
echo "# 3. Annotating protein sequences with Pfam domains        #"
echo "############################################################"
PFAM=${pfam}
PEPTIDES="${CTAT_LIB}ref_annot.pep"

HMM_SCAN="${HMM_EXE} --domtblout '${CTAT_LIB}ref_annot.pfam.dat' --noali --domE 0.001 --seed 42 --qformat fasta --cpu 4 '${PFAM}' '${PEPTIDES}' 1>> '${CTAT_LIB}ref_annot_$(date -I).pfam.log'"
if [ ! -f "${CTAT_LIB}ref_annot.pfam.dat" ] ; then
    # Build pfam domain data
    echo "Running $HMM_SCAN"
    if $HMM_SCAN ; then
        echo "Pfam domain index built"
    else
        echo "Annotation failed; see error message"
    fi
else
    echo -e "${CTAT_LIB}ref_annot.pfam.dat has already been generated..."
fi
