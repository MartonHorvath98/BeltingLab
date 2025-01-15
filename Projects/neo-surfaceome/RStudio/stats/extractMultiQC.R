################################################################################
# 1.) Set up the work environment and the directory structure                  #
################################################################################
# Set downstream path
## folders for: data and figures
data_dir <- "data" # general for every results (only created once)
data_needed <- FALSE # flag to check if the data human.folder is needed
if (!dir.exists(file.path(data_dir))) {
  data_needed <- TRUE
  dir.create(file.path(data_dir)) # create the data human.folder
}

# Figures' folder recreated with current date
date <- format(Sys.Date(), "%Y-%m-%d") # get the current date
if (!dir.exists(file.path(date))) {
  dir.create(file.path(date)) # create the dated results folder
}

# Figure directory
plots <- "plots"
if (!dir.exists(file.path(date, plots))) {
  dir.create(file.path(date, plots)) # create the plots folder
}

################################################################################
# 2.) Load the data and the required files                                     #
################################################################################
# Load the user-defined functions
source("./packages.R")
source("./functions.R")
cat(crayon::white$bold("Loading the data and the required files\n"))

if (!file.exists(
  paste(file.path(human.folder, data_dir), "total_readcounts.xlsx", sep = "/")
)) {
  # EXTRACT ANNOTATIONS
  cat(crayon::white("-> Extracting sections from mutliqc_data.json...\n"))
  # Extracting transcript annotations from the GFF3 file into a TxDb object
  human.txdb <- makeTxDbFromGFF("../reference/Homo_sapiens.GRCh38.96.gff3", 
                                dataSource = "Ensembl", 
                                organism = "Homo sapiens")
  # Extracting gene features from the TxDb object
  human.genes <- exonsBy(human.txdb, by = "gene")
  
  # READ IN BAM ALIGNMENT FILES
  cat(crayon::white("-> Extracting read counts from the BAM files...\n"))
  human.path  <- choose.dir(getwd(), "Select the directory containing the '.bam' files")
  human.files <- list.files(human.path, pattern = ".bam$", full.names = T)
  human.list <- BamFileList(human.files, yieldSize = 2000000)
  
  # CALCULATE READCOUNTS
  cat(crayon::white("-> Calculating read counts..."))
  # Set up parallel computing
  register(SnowParam())
  # Calculate read counts
  human.reads <- summarizeOverlaps(features = human.genes, 
                                   reads = human.list, 
                                   # (default) count reads overlapping single exon
                                   mode = "Union",  
                                   # use strand-specificity
                                   ignore.strand = F, 
                                   # paired-end reads
                                   singleEnd = F, 
                                   fragments = T, 
                                   # strand-specificity => first strand, RF
                                   preprocess.reads = invertStrand)
  
  # CREATE DATA FRAME WITH READ COUNTS
  cat(crayon::white("-> creating data frame with readcounts...\n"))
  # Extract count matrix
  human.counts <- assay(human.reads)
  human.counts <- human.counts %>%
    as.data.frame() %>%
    # rename columns
    setNames(str_remove(colnames(.), ".bam"))
  # Save the read counts to an excel file
  write.xlsx(human.counts, 
             paste(file.path(human.folder, data_dir),"total_readcounts.xlsx", sep = "/"),
             keepNA = T, na.string = "NA", colNames = T, rowNames = T)
}