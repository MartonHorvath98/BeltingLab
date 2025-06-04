#!/bin/bash
# Using getopt

######################################################
# Marton Horvath, APR 2025 - last modified: MAJ 2025 #
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

# DOWNLOAD THE UNIPROT TAB FILE
echo "#############################################"
echo "# 1. Download ID mapping tab file (Uniprot) #"
echo "#############################################"
# Set up the URL to Uniprot by organism download link
uniprot_url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/"
# Create directory
output_dir="${output}"
mkdir -p ${output_dir}

# Download the idmapping data files which are updated in conjunction with the UniProt Knowledgebase (UniProtKB).
UNIPROT_ID="HUMAN_9606_idmapping.dat.gz"
if [ ! -f "${output_dir}${UNIPROT_ID}" ] ; then
	get ${string_url}${UNIPROT_ID} ${output_dir} || (echo "Error getting ${UNIPROT_ID}" && exit 1)
	gunzip ${output_dir}${UNIPROT_ID} || (echo "Error unzipping ${UNIPROT_ID}" && exit 1)
    mv ${output_dir}${UNIPROT_ID%.gz} ${output_dir}${UNIPROT_ID%.gz}.txt
else
	echo "${UNIPROT_ID%.gz} have been dowloaded already..."
fi
UNIPROT_ID="${output_dir}${UNIPROT_ID%.gz}.txt"
echo -e "Downloaded ID mapping list: ${UNIPROT_ID}"

# PARSE THE IDMAPPING DATA FILES
echo "#################################"
echo -e "# Processing IDMAPPING files #"
echo "#################################"
# Create the interaction table
input_file="${UNIPROT_ID}"
python3 parse_uniprot.py -i ${input_file} -o ${output_dir}"uniprot_ID_mapping.txt" || (echo "Error processing ${input_file}" && exit 1)

echo -e "Processed interactions saved to: ${output_dir}uniprot_ID_mapping.txt"
rm ${input_file} || (echo "Error removing ${input_file}" && exit 1)
