####################################################
# 1.) Set up working directory and directory tree  #
####################################################
# Set downstream path
## Create the sub folders for: results, data, and pictures
data <- "data"
if (!dir.exists(file.path(data))) {
  dir.create(file.path(data)) # create the data folder
} 
#plots directory
plots <- "plots"
#results directory
results <- "results"
# get the current date
date <- format(Sys.Date(), "%Y-%m-%d")
if (!dir.exists(file.path(date))) {
  dir.create(file.path(date)) # create the dated results folder
  dir.create(file.path(date, results)) # create the results folder
  dir.create(file.path(date, plots)) # create the plots folder
}

#######################################
# 2.) Load functions for the analyses #
#######################################
source("functions.R")
source("packages.R")

############################## 
# 3.) Load analysis data     #
##############################
if (exists("readcounts") == F) {
  mrna.path <- file.path("../../results/06_counts/")
  mrna.files <- list.files(mrna.path, pattern = "counts.txt$", full.names = T)
  file.copy(from = mrna.files, to = file.path(data))
  
  mrna.reads <- lapply(list.files(file.path(data), 
                                  pattern = "counts.txt$",
                                  full.names = T), 
                       # read csv: with header, tab-delimited, skip first row =>
                       # contains quantification parameters
                       read.csv, header = T, sep = "\t", skip = 1)
  
  mrna.counts <- lapply(mrna.reads, 
                        function(x) { x  %>%
                            dplyr::select(c(1,7)) %>%
                            setNames(c("Geneid","Counts")) %>% 
                            dplyr::mutate(Geneid = stringr::str_remove(Geneid, "gene:"))}
  )
  
  mrna.counts <- merge.rec(mrna.counts, by = "Geneid",  all = T, suffixes = c("",""))
  file.names <- list.files(file.path(data),
                           pattern = "counts.txt$", full.names = F)
  
  mrna.counts <- mrna.counts %>%
    # set the gene ID as row names
    tibble::column_to_rownames("Geneid") %>% 
    # clean up the column names
    setNames(str_remove(file.names, ".counts.txt"))
  
  write.table(mrna.counts, paste(data,"readcounts.csv", sep = "/"), 
              sep =",", na = "NA", dec = ".", row.names = T, col.names = T)
}

##############################
# 4.) Ready count tables     #
##############################
readcounts <- read.csv(file = file.path(data,"readcounts.csv"),
                       sep = ",", header = T, na.strings = NA, row.names = 1)
readcounts <- as.matrix(readcounts)

##############################
# 5.) Prepare metadata table #
##############################
samples <- colnames(readcounts)
coldata <- data.frame("samplenames" = samples) %>%
  # extract the cell line from the sample names
  dplyr::mutate(samplenames = as.factor(samplenames),
                patient = dplyr::case_when(
                  stringr::str_detect(samples, "593") ~ "593",
                  stringr::str_detect(samples, "673") ~ "673"
                  ),
                patient = factor(patient, levels = c("593","673"))) %>%
  # extract the treatment from the sample names
  dplyr::mutate(dimension = dplyr::case_when(
                  stringr::str_detect(samples, ".2D") ~ "2D",
                  stringr::str_detect(samples, ".3D") ~ "3D",
                  TRUE ~ "Tumour"),
                dimension = factor(dimension, levels = c("2D","3D","Tumour"))) %>% 
  # extract the replicates from the sample names
  dplyr::mutate(condition = dplyr::case_when(
                  stringr::str_detect(samples, "H.") ~ "hypoxia",
                  stringr::str_detect(samples, "N.") ~ "normoxia",
                  TRUE ~ "physioxia"),
                condition = factor(condition, levels = c("normoxia","physioxia","hypoxia"))) %>%
  # extract the replicates from the sample names
  dplyr::mutate(run = dplyr::case_when(
    stringr::str_detect(samples, ".1$") ~ "run1",
    stringr::str_detect(samples, ".2$") ~ "run2",
    stringr::str_detect(samples, ".3$") ~ "run3"
  ),
  run = factor(run, levels = c("run1","run2","run3"))) %>%
  # create the experiment and full setup
  dplyr::mutate(groups = as.factor(str_glue("{dimension}_{condition}",
                                               cells = cells,
                                               treatment = treatment)))

deseq <- make_deseq(matrix = readcounts,
                    coldata = coldata,
                    design = "groups")

pca <- make_pca(deseq$rld, group = "condition", labs = c("normoxia","physioxia","hypoxia"), 
                     cols = c("steelblue","salmon","darkred"))
ggsave(paste(file.path(date, plots),"pca.png",sep="/"),
       plot = pca, width = 10, height = 10, units = 'in')

# Create a heatmap of the top 100 genes
top100 <- names(rowVars(assay(deseq$rld))[order(rowVars(assay(deseq$rld)), decreasing = T)[1:100]])
