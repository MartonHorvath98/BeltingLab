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

# DOWNLOAD THE STRING FLAT FILE
echo "##################################################"
echo "# 1. Download interacting protein pairs (STRING) #"
echo "##################################################"
# Set up the URL to the current flat file from STRINGdb
string_url="https://stringdb-downloads.org/download/"
# Create directory
output_dir="${output}"
mkdir -p ${output_dir}

# Download the protein network data (full network, incl. subscores per channel) from STRINGdb
STRING_LINKS="protein.links.full.v12.0/9606.protein.links.detailed.v12.0.txt.gz"
if [ ! -f "${output_dir}${STRING_LINKS}" ] ; then
	get ${string_url}${STRING_LINKS} ${output_dir} || (echo "Error getting ${STRING_LINKS}" && exit 1)
	gunzip ${output_dir}${STRING_LINKS} || (echo "Error unzipping ${STRING_LINKS}" && exit 1)
    mv ${output_dir}${STRING_LINKS%.gz} ${output_dir}${STRING_LINKS%.gz}.txt
else
	echo "${STRING_LINKS%.gz} have been dowloaded already..."
fi
STRING_LINKS="${output_dir}${STRING_LINKS%.gz}.txt"
echo -e "Downloaded interacting protein pairs: ${STRING_LINKS}"

# Download the informations of STRING proteins, incl. their display names and descriptions
STRING_INFO="protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
if [ ! -f "${output_dir}${STRING_INFO}" ] ; then
	get ${string_url}${STRING_INFO} ${output_dir} || (echo "Error getting ${STRING_INFO}" && exit 1)
	gunzip ${output_dir}${STRING_INFO} || (echo "Error unzipping ${STRING_INFO}" && exit 1)
    mv ${output_dir}${STRING_INFO%.gz} ${output_dir}${STRING_INFO%.gz}.txt
else
	echo "${STRING_INFO%.gz} have been dowloaded already..."
fi
STRING_INFO="${output_dir}${STRING_INFO%.gz}.txt"
echo -e "Downloaded protein informations: ${STRING_INFO}"

# PARSE THE STRING DATA FILES
echo "#################################"
echo -e "# Processing PPI from STRING-DB #"
echo "#################################"
# Create the interaction table
input_file="${STRING_LINKS}"
python3 parse_string.py -i ${input_file} -o ${output_dir}"protein_links_clean.tsv" 0,1 '\t' || (echo "Error processing ${input_file}" && exit 1)

echo -e "Processed interactions saved to: ${output_dir}protein_links_clean.tsv"
rm ${input_file} || (echo "Error removing ${input_file}" && exit 1)

echo "##########################################"
echo -e "# Processing protein info from STRING-DB #"
echo "##########################################"
# Create the protein info table
input_file="${STRING_INFO}"
python3 parse_string.py -i ${input_file} -o ${output_dir}"protein_info_clean.tsv" 0 '\t' || (echo "Error processing ${input_file}" && exit 1)

echo -e "Processed interactions saved to: ${output_dir}protein_info_clean.tsv"
rm ${input_file} || (echo "Error removing ${input_file}" && exit 1)

