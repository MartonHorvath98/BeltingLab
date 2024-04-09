#! bin/bash
# Marton Horvath, April 2024

# Create environment
mamba env create -f config/environment.yml

#Set up directory tree
mamba activate glioblastoma

# Run the QC pipeline
echo "############################"
echo "Subsampling raw data files #"
echo "############################"


# Raw data directory
raw_data="data/00_raw"
# Sampled data directory
sample_data="data/01_sample"
mkdir -p $sample_data

# Subsample the raw data
source bin/getSamples.sh -i $raw_data -c "config/" -o $sample_data


# Run the QC pipeline
echo "####################################"
echo "Running quality check pipeline      "
echo "####################################"

# Create the QC results directory
qc_foder="results/01_QC"
mkdir -p $qc_folder

# Run the QC pipeline
source bin/QC.sh -d $raw_data -o $qc_folder

# Run the gene prediction pipeline
echo "####################################"
echo "Running the gene prediction pipeline"
echo "####################################"

# Create the results directory
gene_folder="results/02_gene_prediction"
mkdir -p $gene_folder

# minimum contig length of 1000 bp for the first run
source bin/genePred.sh -d $raw_data -o ${gene_folder}/raw -m 1000 

# Clean the Haemoproteus tartakovskyi genome assembly
echo "###############################################"
echo "Cleaning the Haemoproteus tartakovskyi assembly"
echo "###############################################"

# Create the results directory
clean_data="data/02_clean"
uniprot_db="data/01_uniprot/aves_taxid.tsv"
mkdir -p $clean_data

# Run the cleanAves.sh script
source bin/cleanAves.sh -fa ${raw_data}/Haemoproteus_tartakovskyi.genome \
    -g ${gene_folder}/genemark.Ht.gtf -db $uniprot_db \
    -gc 27 -o $clean_data

# Copy the rest of the raw genomes to the clean data folder
cp -l ${raw_data}/$(ls --hide=Haemoproteus_tartakovskyi.genome ${raw_data}/) ${clean_data}

# Ru-run gene prediction with minimum contig length of 3000 bp
source bin/genePred.sh -d $clean_data -o ${gene_folder}/clean -m 3000 

# Run BUSCO check and get orthologous gene groups
echo "###############################################"
echo " Get orthologous gene groups and check BUSCO   "
echo "###############################################"

# Create the results directory
busco_db="data/03_busco_downloads/apicomplexa_odb10"
busco_folder="results/04_BUSCO"
ortho_folder="results/05_ORTHO"

mkdir -p $busco_folder
mkdir -p $ortho_folder

# Run the getOrtho.sh script
source bin/getOrtho.sh -i $clean_data -g $gene_folder/clean -d $busco_db -b $busco_folder -o $ortho_folder

# Generate the phylogenetic tree with clustalo and iqtree
echo "###############################################"
echo " Generate the phylogenetic tree with iqtree   "
echo "###############################################"

# Create the results directory
phylo_folder="results/06_PHYLO"
mkdir -p $phylo_folder

# Run the createPhyloTree.sh script
source bin/createPhylo_BUSCO.sh -i $busco_folder -p $ortho_folder -o $phylo_folder
source bin/createPhylo_proteinortho.sh -i $ortho_folder -p "Malaria_phylo.proteinortho.tsv" -o $phylo_folder

# Create consensus tree with phylip consense
echo "###############################################"
echo " Create consensus tree with phylip consense   "
echo "###############################################"

# Create the results directory
consensus_folder="results/07_CONTREE"
mkdir -p $consensus_folder

# Run the consensusTree.sh script
source bin/consensusTree.sh -i ${phylo_folder}/tree -o ${consensus_folder}/complete -m false
source bin/consensusTree.sh -i ${phylo_folder}/tree_proteinortho -o ${consensus_folder}/proteinortho -m false
source bin/consensusTree.sh -i ${phylo_folder}/noTG_tree -o ${consensus_folder}/minOut -m true

