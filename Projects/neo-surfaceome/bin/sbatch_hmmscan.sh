#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p node -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J domain-map
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
 
# Marton Horvath, MAR 2025 - last modified: MAR 2025
# ------------ Start a screen environment ----- #
screen -S PfamScan
# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools hmmer/3.3.2
# ------------ Set global variables ------------ #
# Path to the PFAM database
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/
PFAM="${WHARF}resources/database/pfam/Pfam-A.hmm"
# ENSEMBL protein sequences fasta
FASTA="${WHARF}resources/database/ensembl/grch38_v103_proteins.fasta"
# Output directory
OUTPUT_DIR="${WHARF}resources/database/"

# ----------- Check executables ---------------- #
## 1. HMMSCAN
HMM_EXE=./hmmscan
if [ ! -x "$HMM_EXE" ] ; then
	if ! which hmmscan; then
		echo "Could not find 'hmmscan' in current directory or in PATH"
		exit 1
	else
		HMM_EXE=`which hmmscan`
	fi
fi

# echo "##############################################"
# echo "# Splitting the fasta file per protein  #"
# echo "##############################################"

# export TMP_DIR="${OUTPUT_DIR}ensembl/tmp/"

# if [ ! -d ${TMP_DIR} ] ; then 
#   mkdir -p ${TMP_DIR}
    
    # Create ENSEMBL protein FASTA files
#    perl -ane ' if(/^>(ENSP\d+)/) { $file="$ENV{TMP_DIR}$1.fasta"; if (-e $file) { warn "Skipping existing file: $file\n"; next } close O if $.>1; open(O, ">", $file) or die "Cannot write to $file: $!"; print O $_ } else { print O $_ }' ${FASTA}
# fi

echo "############################################################"
echo "# Annotating ENSEMBL's protein sequences with Pfam domains #"
echo "############################################################"

# if [ -d ${TMP_DIR} ] ; then
#    for seq in ${TMP_DIR}*.fasta; do
srun --exclusive --ntasks=4 --cpus-per-task=4 "${HMM_EXE}" --domtblout "${OUTPUT_DIR}pfam/grch38_prot_pfam.dat" \
    --noali --domE 0.001 --seed 42 --qformat fasta --cpu 4 ${PFAM} ${FASTA} 1>> "${OUTPUT_DIR}pfam/grch38_prot_pfam_$(date -I).log"
# done
# fi
# rm -r ${TMP_DIR]

# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
