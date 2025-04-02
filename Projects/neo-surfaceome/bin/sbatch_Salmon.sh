#!/bin/bash -l 
#SBATCH -A sens2020018
#SBATCH -p node -n 1
#SBATCH -t 3-00:00:00
#SBATCH -J align-Salmon
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4

# Marton Horvath, Mar 2025; last modified: Mar 2025
# ------------ Start a screen environment ----- #
screen -S SalmonQuants
# ------------ Load needed modules ------------ #
# Load bioinfo-tools with programs
module load bioinfo-tools Salmon/1.10.1 samtools/1.20

# ------------ Export paths -------------------- #
export WHARF=/proj/sens2020018/nobackup/wharf/hmarton/hmarton-sens2020018/
export CTAT=${WHARF}refs/ctat/

# ------------ Set global variables ------------ #
SAMPLES="$(cat config/samples.csv | cut -d ',' -f1 | sed 's/Sample_//g')"
INPUT_DIR="${WHARF}results/02_trim/"
eval echo "STAR alignments from: ${INPUT_DIR}"

REF="${CTAT}ref_genome.fa.salmon.idx/"
GENOME="${CTAT}ref_genome.fa"
TRANSCRIPTOME="${CTAT}ref_annot.cds"
GTF="${CTAT}ref_annot.gtf"

SALMON_DIR="${WHARF}results/04_counts/"

# ------- Check command executable paths -------
## 1. Salmon
SALMON_EXE=./salmon
if [ ! -x "$SALMON_EXE" ] ; then
    if ! which salmon ; then
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

# ------- Build Salmon index (if missing) -------
if [ ! -d "${REF}" ] ; then
    echo "##########################################################"
    echo "# Build Salmon index                                     #"
    echo "##########################################################"
    # Create decoys.txt file
    mkdir -p "${REF}"
    grep "^>" ${GENOME} | cut -d " " -f 1 > "${CTAT}decoys.txt"
    sed -i.bak -e 's/>//g' "${CTAT}decoys.txt"
    
    # Concatenate transcriptome and genome
    cat ${TRANSCRIPTOME} ${GENOME} > "${CTAT}gentrome.fa"

    # Build Salmon index
    SALMON_INDEX="srun --exclusive --ntasks=1 --cpus-per-task=4 ${SALMON_EXE} index -t ${CTAT}gentrome.fa -d ${CTAT}decoys.txt -i ${REF} -k 31 --gencode"
    echo Running $SALMON_INDEX
    if $SALMON_INDEX ; then
	# Remove concatenated transcriptome and genome
	rm -f "${CTAT}decoys.txt.bak" "${CTAT}gentrome.fa"
	echo "Salmon index built"
    else
	echo "Index building failed; see error message"
    fi
fi

# ------- Run Salmon Quantification ------- #
TMP_DIR="${SALMON_DIR}tmp/"
mkdir -p "${TMP_DIR}"

for sample in ${SAMPLES}
do
    # Read the trimmed forward and reverse reads
    #echo "Extract aligned features from STAR bam: ${sample}"
    #srun --exclusive --ntasks=1 --cpus-per-task=8 ${SAMTOOLS_EXE} fastq -@ 8 \
	#${INPUT_DIR}${sample}-Aligned.sortedByCoord.out.bam \
	#-1 ${TMP_DIR}${sample}_R1.fq.gz \
	#-2 ${TMP_DIR}${sample}_R2.fq.gz \
	#-0 /dev/null -s /dev/null -n

    # Quantify STAR alignments with Salmon
    srun --exclusive --ntasks=2 --cpus-per-task=4 $SALMON_EXE quant -i "${REF}" -l "ISR" \
	-1 "${INPUT_DIR}${sample}_R1_001_val_1.fq.gz" \
        -2 "${INPUT_DIR}${sample}_R2_001_val_2.fq.gz" \
        -p 4 -g "${GTF}" --seqBias -o "${SALMON_DIR}${sample}"

    #rm -r ${TMP_DIR}
done

# ------------ Deactivate the environment ------------ #
# Workflow finished
echo "####################################"
echo " Workflow finished                  "
echo "####################################"

exit
