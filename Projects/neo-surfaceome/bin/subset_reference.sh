#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p core -n 4
#SBATCH -t 1-00:00:00
#SBATCH -J subset-reference

# Marton Horvath, DEC 2024 - last modified: DEC 2024
# ------------ Start a screen environment ----- #
screen -S SubsetRefGenome
# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools star/2.7.8a BEDTools/2.31.1 samtools/1.20

# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/

INPUT_DIR=${WHARF}refs/ctat/
OUTPUT_DIR=${WHARF}refs/ctat_sub/
ALIGN_DIR=${WHARF}results/03_align/
SNP_DIR=${WHARF}results/05_SNP/

# ------------ Set global variables ------------ #
# Set the number of cores
CORES=4
eval echo "Number of cores: ${CORES}"# Set the number of cores
# Export sample names
SAMPLES="$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')"
# Export gene list
GENES=${WHARF}config/surfme_genes.txt

# ------- Check command executable paths -------
## 2. Samtools for indexing BAM files
SAMTOOLS_EXE=./samtools
if [ ! -x "$SAMTOOLS_EXE" ] ; then
	if ! which samtools ; then
		echo "Could not find samtools in current directory or in PATH"
		exit 1
	else
		SAMTOOLS_EXE=`which samtools`
	fi
fi

REF="ref_genome.fa"
GTF="ref_annot.gtf"
BED="genes.bed"

echo "###############################################"
echo "# Subset reference genome                     #"
echo "###############################################"
# ------------ Extract sequence and annotation for surface genes ---- #
# Filter GTF file
grep -Ff ${GENES} ${INPUT_DIR}${GTF} > ${OUTPUT_DIR}${GTF}

#Convert GTF to BED
awk '$3 == "gene"' ${OUTPUT_DIR}${GTF} | awk '{ print $1, $4-1, $5, $10 }' OFS="\t" > ${OUTPUT_DIR}${BED}

echo "###############################################"
echo "# Subset BAM files                            #"
echo "###############################################"
# ------------ Extract features defined by the BED file ---- #
SAMTOOLS_FILTER="${SAMTOOLS_EXE} view -b -h -L"

if [ -d "$ALIGN_DIR" ] ; then
    # Iterate each sample
    for sample in ${SAMPLES}
    do
	echo -e "Subsetting: ${sample}"
	${SAMTOOLS_FILTER} ${OUTPUT_DIR}${BED} ${ALIGN_DIR}${sample}-Aligned.sortedByCoord.out.bam --threads ${CORES} > ${SNP_DIR}${sample}-SURFME.sortedByCoord.out.bam
	
    done
fi
       
# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit

