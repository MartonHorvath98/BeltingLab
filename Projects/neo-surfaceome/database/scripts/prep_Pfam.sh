#!/bin/bash
# Using getopt

#####################################################
## Marton Horvath, MAR 2025 - last modified: APR 2025
#####################################################

# ------------ Custom functions ------------ #
# Usage function
usage() {
    echo "Usage: $0 -o <output_dir [PATH]>"
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
while getopts "o:" opt; do
    case "$opt" in
		o) output=$OPTARG ;;
        *) usage ;;
    esac
done

# Create and activate the environment
mamba env create -f config/buildenv.yml

# ------- Check command executable paths ------ #
HMMER_EXE=./hmmpress
if [ ! -x "$HMMER_EXE" ] ; then
	if ! which hmmpress ; then
		echo "Could not find hmmer3 in current directory or in PATH"
		exit 1
	else
		HMMER_EXE=`which hmmpress`
	fi
fi

# DOWNLOAD THE PFAM HMM LIBRARY AND ACTIVE SITE TARBALL
echo "#################################################"
echo "# 1. Download Pfam-A.hmm and active_site.dat... #"
echo "#################################################"
# Set up the URL to the current release of Pfam-A.hmm
pfam_url="https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/"
# Create directory structure for Ensembl protein FASTA files
output_dir="${output}"
mkdir -p ${output_dir}

# Download the Pfam HMM library for Pfam-A families (unless it already exists)
PFAM_A="Pfam-A.hmm.gz"
if [ ! -f "${output_dir}${PFAM_A}" ] ; then
	get ${pfam_url}${PFAM_A} ${output_dir} || (echo "Error getting ${PFAM_A}" && exit 1)
	gunzip ${output_dir}${PFAM_A} || (echo "Error unzipping ${PFAM_A}" && exit 1)
else
	echo "${PFAM_A%.gz} have been dowloaded already..."
fi
PFAM_A="${output_dir}${PFAM_A%.gz}"
echo -e "Downloaded Pfam HMM library for Pfam-A families: ${PFAM_A}"

# Download the active_site.dat file (unless it already exists)
ACTIVE_SITE="active_site.dat.gz"
if [ ! -f "${output_dir}${ACTIVE_SITE}" ] ; then
	get ${pfam_url}${ACTIVE_SITE} ${output_dir} || (echo "Error getting ${ACTIVE_SITE}" && exit 1)
	gunzip ${output_dir}${ACTIVE_SITE} || (echo "Error unzipping ${ACTIVE_SITE}" && exit 1)
else
	echo "${ACTIVE_SITE%.gz} have been dowloaded already..."
fi
ACTIVE_SITE="${output_dir}${ACTIVE_SITE%.gz}"
echo -e "Downloaded active_site.dat file: ${ACTIVE_SITE}"

# generate binary files for Pfam-A.hmm (Pfam-A.hmm, Pfam-A.hmm.dat, active_site.dat files need to be downloaded)
echo "##############################################"
echo "# 2. Generating binary index for Pfam-A.hmm  #"
echo "##############################################"
# Create binary index
extensions=("*.h3f" "*.h3i" "*.h3m" "*.h3p")

for ext in "${extensions[@]}"; do
    if ! ls ${output_dir}$ext &>/dev/null; then
        echo "Missing file with extension: $ext"
		HMMER_INDEX="${HMMER_EXE} ${output_dir}PFam-A.hmm"
		if "${HMMER_INDEX}" ; then
			echo "Pfam index built"
		else
			echo "Index building failed; see error message"
		fi
	fi
done
echo -e "Pfam index have been generated already..."