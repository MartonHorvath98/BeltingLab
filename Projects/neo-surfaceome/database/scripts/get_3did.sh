#!/bin/bash
# Using getopt

######################################################
# Marton Horvath, MAR 2025 - last modified: APR 2025 #
######################################################

# ------------ Custom functions ------------ #
# Usage function
usage() {
    echo "Usage: $0 -i <input_file>"
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

# DOWNLOAD THE 3DID FLAT FILE
echo "###############################################"
echo "# 1. Download interacting domain pairs (3DID) #"
echo "###############################################"
# Set up the URL to the current flat file of 3DID
did_url="https://3did.irbbarcelona.org/download/current/"
# Create directory structure for Ensembl protein FASTA files
output_dir="${output}"
mkdir -p ${output_dir}

# Download the flat file from 3DID containing interacting domain pairs (ID) 
FLAT_3DID="3did_flat.gz"
if [ ! -f "${output_dir}${FLAT_3DID}" ] ; then
	get ${did_url}${FLAT_3DID} ${output_dir} || (echo "Error getting ${FLAT_3DID}" && exit 1)
	gunzip ${output_dir}${FLAT_3DID} || (echo "Error unzipping ${FLAT_3DID}" && exit 1)
    mv ${output_dir}${FLAT_3DID%.gz} ${output_dir}${FLAT_3DID%.gz}.dat
else
	echo "${FLAT_3DID%.gz} have been dowloaded already..."
fi
FLAT_3DID="${output_dir}${FLAT_3DID%.gz}.dat"
echo -e "Downloaded interacting domain pairs: ${FLAT_3DID}"

# PARSE THE 3DID FLAT FILE
echo "###################################################"
echo -e "# Processing Pfam interactions from 3did_flat.dat #"
echo "###################################################"
# Create the interaction table
input_file="${FLAT_3DID}"
preprocessed_file="domain_interactions.tsv"
less "$input_file" | grep "^#=ID" | cut -f4,5 | perl -ane '
	BEGIN { print "domain_A\tdomain_B\n" }
    $F[0] =~ s/.*(PF\d+).*/$1/;
    $F[1] =~ s/.*(PF\d+).*/$1/;
    print "$F[0]\t$F[1]\n$F[1]\t$F[0]\n";
' | sort -u > "${output_dir}/${preprocessed_file}"

echo -e "Processed interactions saved to: ${output_dir}/${preprocessed_file}"
rm ${FLAT_3DID} || (echo "Error removing ${FLAT_3DID}" && exit 1)
