#!/bin/bash
# Using getopt

#####################################################
## Marton Horvath, JUN 2024 - last modified: MAR 2025
#####################################################
# Usage function
usage() {
    echo "Usage: $0 -r <ensembl_release [INT]> -p <pfam_lib [PATH]> -f <fasta [STR] -o <output_dir [PATH]>"
    exit 1
}
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

# Parse command-line options
while getopts "r:p:f:o:" opt; do
    case "$opt" in
        r) ensembl_release=$OPTARG ;;
		p) pfam_lib=$OPTARG ;;
		f) fasta=$OPTARG ;;
		o) output=$OPTARG ;;
        *) usage ;;
    esac
done

# QC for user arguments
if [[ -z "$ensembl_release" || ! "$ensembl_release" =~ ^[0-9]+$ ]]; then
    echo "Error: Invalid or missing Ensembl release (-r). Must be an integer."
    usage
fi

if [[ -z "$pfam_lib" || ! -d "$pfam_lib" ]]; then
    echo "Error: Invalid or missing Pfam library path (-p). Directory does not exist."
    usage
fi

if [[ -z "$fasta" || ! "$fasta" =~ ^[a-zA-Z0-9._-]+$ ]]; then
    echo "Error: Invalid or missing FASTA filename (-f). Must be a valid string."
    usage
fi


# Create and activate the environment
# mamba env create -f environment.yml && conda activate pfamscan

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

# DOWNLOAD REFERENCE GENOME AND ANNOTATION
echo "##############################################"
echo "# 1. Download Ensembl protein FASTA files... #"
echo "##############################################"
# Set up the Ensembl release version
release=$ensembl_release
# Set up the Ensembl base URLs
ENSEMBL_PROTEIN=ftp://ftp.ensembl.org/pub/release-${release}/fasta/homo_sapiens/pep
# Create directory structure for Ensembl protein FASTA files
output_dir="${output}ensembl_v${release}/"
mkdir -p ${output_dir}

# Download the GRCh38 reference genome fasta file (unless it already exists)
PEPTIDES=Homo_sapiens.GRCh38.pep.all.fa.gz
if [ ! -f "${output_dir}${fasta}" ] ; then
	get ${ENSEMBL_PROTEIN}/$PEPTIDES ${output_dir} || (echo "Error getting $PEPTIDES" && exit 1)
	gunzip ${output_dir}$PEPTIDES || (echo "Error unzipping $PEPTIDES" && exit 1)
	mv ${output_dir}${PEPTIDES%.gz} "${output_dir}${fasta}"
else
	echo "Sequences have been dowloaded already..."
fi

PEPTIDES="${output_dir}${fasta}"
echo -e "Downloaded Ensembl protein sequences are in: ${PEPTIDES}"

# generate binary files for Pfam-A.hmm (Pfam-A.hmm, Pfam-A.hmm.dat, active_site.dat files need to be downloaded)
echo "##############################################"
echo "# 2. Generating binary index for Pfam-A.hmm  #"
echo "##############################################"
# Create STAR index directory
extensions=("*.h3f" "*.h3i" "*.h3m" "*.h3p")

for ext in "${extensions[@]}"; do
    if ! ls ${pfam_lib}$ext &>/dev/null; then
        echo "Missing file with extension: $ext"
		HMMER_INDEX="${HMMER_EXE} ${pfam_lib}PFam-A.hmm"
		if "${HMMER_INDEX}" ; then
			echo "Pfam index built"
		else
			echo "Index building failed; see error message"
		fi
	fi
done
echo -e "Pfam index have been generated already..."

echo "##############################################"
echo "# 3. Generating binary index for Pfam-A.hmm  #"
echo "##############################################"

tmp_dir="${output_dir}tmp/"
mkdir -p ${tmp_dir}
# Create ENSEMBL protein FASTA files
less ${PEPTIDES} | perl -ne ' if(/^>(ENSP\d+)/) {
	$file="resources/ensembl_v113/tmp/$1.fasta"; 
    if (-e $file) { warn "Skipping existing file: $file\n"; next } 
    close O if $.>1; 
    open(O, ">", $file) or die "Cannot write to $file: $!"; 
    print O $_ 
} else { print O $_ }'

# Run pfam_scan.pl. HMM files of Pfam-A have been downloaded to
for s in ${tmp_dir}*.fasta; do perl scripts/pfam_scan.pl -fasta $s -dir ${pfam_lib}/ | tee -a ${output_dir}/ensp_pfam_v${release}.tsv; done
rm -r ${tmp_dir}