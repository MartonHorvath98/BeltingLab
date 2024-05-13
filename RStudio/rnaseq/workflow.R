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
plots_dir <- "plots"
#results directory
results_dir <- "results"
# get the current date
date <- format(Sys.Date(), "%Y-%m-%d")
if (!dir.exists(file.path(date))) {
  dir.create(file.path(date)) # create the dated results folder
  dir.create(file.path(date, results_dir)) # create the results folder
  dir.create(file.path(date, plots_dir)) # create the plots folder
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
# 4.) Load count read counts #
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
                dimension = factor(dimension, levels = c("Tumour","2D","3D"))) %>% 
  # extract the replicates from the sample names
  dplyr::mutate(oxygen = dplyr::case_when(
                  stringr::str_detect(samples, "H.") ~ "hypoxia",
                  stringr::str_detect(samples, "N.") ~ "normoxia",
                  TRUE ~ "physioxia"),
                oxygen = factor(oxygen, levels = c("physioxia","normoxia","hypoxia"))) %>%
  # extract the replicates from the sample names
  dplyr::mutate(run = dplyr::case_when(
    stringr::str_detect(samples, ".1$") ~ "run1",
    stringr::str_detect(samples, ".2$") ~ "run2",
    stringr::str_detect(samples, ".3$") ~ "run3"
  ),
  run = factor(run, levels = c("run1","run2","run3"))) %>%
  # create the experiment and full setup
  dplyr::mutate(
    condition = factor(str_glue("{oxygen}.{dimension}"),
                       levels = c("physioxia.Tumour",
                                  "normoxia.2D","hypoxia.2D", "physioxia.3D")),
    groups = as.factor(str_glue("{patient}_{dimension}_{oxygen}",
                                               cells = cells,
                                               treatment = treatment)))


##############################
# 6.) Run DESeq2 analysis    #
##############################
deseq <- make_deseq(matrix = readcounts,
                    coldata = coldata,
                    design = "condition + condition:patient")

resultsNames(deseq$dds)

invitro.deseq <- make_deseq(matrix = readcounts[,-c(10,20)],
                            coldata = coldata[-c(10,20),] %>% droplevels(.),
                            design = "condition + condition:patient")

# ------------------ Total results ---------------------------------------------
# Extract significant results for the Tumor samples vs 2D normoxia
results.2DN_vs_Tumor <- get_results(dds = deseq$dds, sig_log2FC = 1.5, sig_pval = 0.05,
                             contrast=list("condition_normoxia.2D_vs_physioxia.Tumour"),
                             name = "condition_normoxia.2D_vs_physioxia.Tumour")

# Extract significant results for the Tumor samples vs 2D hypoxia
results.2DH_vs_Tumor <- get_results(dds = deseq$dds, sig_log2FC = 1.5, sig_pval = 0.05,
                             contrast=list("condition_hypoxia.2D_vs_physioxia.Tumour"),
                             name = "condition_hypoxia.2D_vs_physioxia.Tumour")

# Extract significant results for the Tumor samples vs 3D organoids
results.3D_vs_Tumor <- get_results(dds = deseq$dds, sig_log2FC = 1.5, sig_pval = 0.05,
                             contrast=list("condition_physioxia.3D_vs_physioxia.Tumour"),
                             name = "condition_physioxia.3D_vs_physioxia.Tumour")

# Save results to excel files
sapply(c("df","sig_df"), function(x){
  # Tumor samples vs 2D normoxia
  write.csv(results.2DN_vs_Tumor[[x]],
            file.path(date, results_dir, paste0("Tumor_vs_2DN_",x,".csv")),
            row.names = T)
  # Tumor samples vs 2D hypoxia
  write.csv(results.2DH_vs_Tumor[[x]], 
            file.path(date, results_dir, paste0("Tumor_vs_2DH_",x,".csv")),
            row.names = T)
  # Tumour samples vs 3D organoids
  write.csv(results.3D_vs_Tumor[[x]], 
            file.path(date, results_dir, paste0("Tumor_vs_3D_",x,".csv")),
            row.names = T)
})

# --------------------- In-vitro results ---------------------------------------
# Extract significant results for the 2D hypoxia vs 2D normoxia
results.2DH_vs_2DN <- get_results(dds = invitro.deseq$dds, sig_log2FC = 1.5, sig_pval = 0.05,
                             contrast=list("condition_hypoxia.2D_vs_normoxia.2D"),
                             name = "condition_hypoxia.2D_vs_normoxia.2D")

# Extract significant results for the 3D organoids vs 2D normoxia
results.3D_vs_2DN <- get_results(dds = invitro.deseq$dds, sig_log2FC = 1.5, sig_pval = 0.05,
                             contrast=list("condition_physioxia.3D_vs_normoxia.2D"),
                             name = "condition_physioxia.3D_vs_normoxia.2D")


# Save results to excel files
sapply(c("df","sig_df"), function(x){
  # 2D hypoxia vs 2D normoxia
  write.csv(results.2DH_vs_2DN[[x]],
            file.path(date, results_dir, paste0("2DH_vs_2DN_",x,".csv")),
            row.names = T)
  # 3D organoids vs 2D normoxia
  write.csv(results.3D_vs_2DN[[x]], 
            file.path(date, results_dir, paste0("3D_vs_2DN_",x,".csv")),
            row.names = T)
})

################################
# 7.) Overrepresentation (ORA) #
################################
# ---------------------- Total results -----------------------------------------
# 1 KEGG pathways - ORA analysis
## 1.1 Get entrez IDs: gene set of interest and background
genes.2DN_vs_Tumor <- get_genelist(results.2DN_vs_Tumor) # 2D normoxia vs Tumour
genes.2DH_vs_Tumor <- get_genelist(results.2DH_vs_Tumor) # 2D hypoxia vs Tumour
genes.3D_vs_Tumor <- get_genelist(results.3D_vs_Tumor) # 3D organoids vs Tumour

## 1.2 Run the ORA analysis
KEGG.2DN_vs_Tumor <- list()
KEGG.2DN_vs_Tumor$df <- kegg_results(genes.2DN_vs_Tumor$interest,
                                     genes.2DN_vs_Tumor$background)

KEGG.2DH_vs_Tumor <- list()
KEGG.2DH_vs_Tumor$df <- kegg_results(genes.2DH_vs_Tumor$interest,
                                     genes.2DH_vs_Tumor$background)

KEGG.3D_vs_Tumor <- list()
KEGG.3D_vs_Tumor$df <- kegg_results(genes.3D_vs_Tumor$interest,
                                     genes.3D_vs_Tumor$background)

## 1.3 Save the results
openxlsx::write.xlsx(KEGG.2DN_vs_Tumor$df, 
                     file.path(date, results_dir, "2DN_vs_Tumor_KEGG.xlsx"))

openxlsx::write.xlsx(KEGG.2DH_vs_Tumor$df,
                     file.path(date, results_dir, "2DH_vs_Tumor_KEGG.xlsx"))

openxlsx::write.xlsx(KEGG.3D_vs_Tumor$df,
                     file.path(date, results_dir, "3D_vs_Tumor_KEGG.xlsx"))

# -------------------------- In-vitro results -----------------------------------
# 2 KEGG pathways - ORA analysis
## 2.1 Get entrez IDs: gene set of interest and background
genes.2DH_vs_2DN <- get_genelist(results.2DH_vs_2DN) # 2D hypoxia vs 2D normoxia
genes.3D_vs_2DN <- get_genelist(results.3D_vs_2DN) # 3D organoids vs 2D normoxia

## 2.2 Run the ORA analysis
KEGG.2DH_vs_2DN <- list()
KEGG.2DH_vs_2DN$df <- kegg_results(genes.2DH_vs_2DN$interest,
                                   genes.2DH_vs_2DN$background)

KEGG.3D_vs_2DN <- list()
KEGG.3D_vs_2DN$df <- kegg_results(genes.3D_vs_2DN$interest,
                                   genes.3D_vs_2DN$background)

## 2.3 Save the results
openxlsx::write.xlsx(KEGG.2DH_vs_2DN$df, 
                     file.path(date, results_dir, "2DH_vs_2DN_KEGG.xlsx"))
openxlsx::write.xlsx(KEGG.3D_vs_2DN$df, 
                     file.path(date, results_dir, "3D_vs_2DN_KEGG.xlsx"))

####################################################
# 8.) Gene Set Enrichment Analysis (GSEA)          #
####################################################
# ---------------------- Total results -----------------------------------------
## 1.1 GSEA GO analysis
GO.2DN_vs_Tumor <- list()
GO.2DN_vs_Tumor$df <- go_results(genes.2DN_vs_Tumor$interest,
                                 genes.2DN_vs_Tumor$background,
                                 type = 'ALL')
GO.2DH_vs_Tumor <- list()
GO.2DH_vs_Tumor$df <- go_results(genes.2DH_vs_Tumor$interest,
                                 genes.2DH_vs_Tumor$background,
                                 type = 'ALL')

GO.3D_vs_Tumor <- list()
GO.3D_vs_Tumor$df <- go_results(genes.3D_vs_Tumor$interest,
                                 genes.3D_vs_Tumor$background,
                                 type = 'ALL')

# Save the results
openxlsx::write.xlsx(GO.2DN_vs_Tumor$df,
                     file.path(date, results_dir, "2DN_vs_Tumor_GO.xlsx"))

openxlsx::write.xlsx(GO.2DH_vs_Tumor$df,
                     file.path(date, results_dir, "2DH_vs_Tumor_GO.xlsx"))

openxlsx::write.xlsx(GO.3D_vs_Tumor$df,
                     file.path(date, results_dir, "3D_vs_Tumor_GO.xlsx"))

# -------------------------- In-vitro results ----------------------------------
## 1.2 GSEA GO analysis
GO.2DH_vs_2DN <- list()
GO.2DH_vs_2DN$df <- go_results(genes.2DH_vs_2DN$interest,
                               genes.2DH_vs_2DN$background,
                               type = 'ALL')

GO.3D_vs_2DN <- list()
GO.3D_vs_2DN$df <- go_results(genes.3D_vs_2DN$interest,
                               genes.3D_vs_2DN$background,
                               type = 'ALL')

# Save the results
openxlsx::write.xlsx(GO.2DH_vs_2DN$df,
                     file.path(date, results_dir, "2DH_vs_2DN_GO.xlsx"))
openxlsx::write.xlsx(GO.3D_vs_2DN$df,
                     file.path(date, results_dir, "3D_vs_2DN_GO.xlsx"))
####################################################
# 9.) Semantic similarity analysis                 #
####################################################
source("helpers.R")
# ---------------------- Total results -----------------------------------------
ontology <- c("BP","CC","MF")

# 1.1 Prepare gene sets for semantic similarity analysis
SIMMAT.2DN_vs_Tumor <- list()
SIMMAT.2DN_vs_Tumor <- GO.2DN_vs_Tumor$df %>%
  dplyr::group_split(ONTOLOGY)
SIMMAT.2DN_vs_Tumor <- setNames(SIMMAT.2DN_vs_Tumor, ontology)

SIMMAT.2DH_vs_Tumor <- list()
SIMMAT.2DH_vs_Tumor <- GO.2DH_vs_Tumor$df %>%
  dplyr::group_split(ONTOLOGY)
SIMMAT.2DH_vs_Tumor <- setNames(SIMMAT.2DH_vs_Tumor, ontology)

SIMMAT.3D_vs_Tumor <- list()
SIMMAT.3D_vs_Tumor <- GO.3D_vs_Tumor$df %>%
  dplyr::group_split(ONTOLOGY)
SIMMAT.3D_vs_Tumor <- setNames(SIMMAT.3D_vs_Tumor, ontology)

# 1.2 Run the semantic similarity analysis
## 2D normoxia vs Tumour
SIM.2DN_vs_Tumor <- make_simMatrix(SIMMAT.2DN_vs_Tumor) # make similarity matrix
SCORES.2DN_vs_Tumor <- lapply(SIMMAT.2DN_vs_Tumor, function(x) setNames(-log10(x$pvalue), x$ID))
for(var in names(SCORES.2DN_vs_Tumor)){
  SCORES.2DN_vs_Tumor[[var]] <- SCORES.2DN_vs_Tumor[[var]][which(names(SCORES.2DN_vs_Tumor[[var]]) %in% row.names(SIM.2DN_vs_Tumor[[var]]))]
}  # calculate the scores for the significant terms
REDUCED.2DN_vs_Tumor <- list()
for(var in names(SIM.2DH_vs_Tumor)){
  REDUCED.2DN_vs_Tumor[[var]] <- get_reducedTerms(simm = SIM.2DN_vs_Tumor[[var]], 
                                         scores = SCORES.2DN_vs_Tumor[[var]], 
                                         0.9)
} # reduce collinearity

## 2D hypoxia vs Tumour
SIM.2DH_vs_Tumor <- make_simMatrix(SIMMAT.2DH_vs_Tumor) # make similarity matrix
SCORES.2DH_vs_Tumor <- lapply(SIMMAT.2DH_vs_Tumor, function(x) setNames(-log10(x$pvalue), x$ID))
for (var in names(SCORES.2DH_vs_Tumor)){
  SCORES.2DH_vs_Tumor[[var]] <- SCORES.2DH_vs_Tumor[[var]][which(names(SCORES.2DH_vs_Tumor[[var]]) %in% row.names(SIM.2DH_vs_Tumor[[var]]))]
} # calculate the scores for the significant terms
REDUCED.2DH_vs_Tumor <- list()
for(var in names(SIM.2DH_vs_Tumor)){
  REDUCED.2DH_vs_Tumor[[var]] <- get_reducedTerms(simm = SIM.2DH_vs_Tumor[[var]], 
                                         scores = SCORES.2DH_vs_Tumor[[var]], 
                                         0.9)
} # reduce collinearity

## 3D organoids vs Tumour
SIM.3D_vs_Tumor <- make_simMatrix(SIMMAT.3D_vs_Tumor) # make similarity matrix
SCORES.3D_vs_Tumor <- lapply(SIMMAT.3D_vs_Tumor, function(x) setNames(-log10(x$pvalue), x$ID))
for(var in names(SCORES.3D_vs_Tumor)){
  SCORES.3D_vs_Tumor[[var]] <- SCORES.3D_vs_Tumor[[var]][which(names(SCORES.3D_vs_Tumor[[var]]) %in% row.names(SIM.3D_vs_Tumor[[var]]) )]
} # calculate the scores for the significant terms
REDUCED.3D_vs_Tumor <- list()
for(var in names(SIM.3D_vs_Tumor)){
  REDUCED.3D_vs_Tumor[[var]] <- get_reducedTerms(simm = SIM.3D_vs_Tumor[[var]], 
                                         scores = SCORES.3D_vs_Tumor[[var]], 
                                         0.9)
} # reduce collinearity

## 1.3 Create compound data frame for visualization
## 2D normoxia vs Tumor
COMP.2DN_vs_Tumor <- bind_df(SIMMAT.2DN_vs_Tumor, KEGG.2DN_vs_Tumor$df)
COMP.2DN_vs_Tumor$category <- factor(COMP.2DN_vs_Tumor$category, 
                                     levels=c("BP", "MF", "CC", "KEGG"))

COMP.2DN_vs_Tumor <- make_GObase(COMP.2DN_vs_Tumor, results.2DN_vs_Tumor$sig_df)
COMP.2DN_vs_Tumor$slim <- ifelse(COMP.2DN_vs_Tumor$id %in% 
                                   c(processes,
                                     REDUCED.2DN_vs_Tumor$MF$subset$go,
                                     REDUCED.2DN_vs_Tumor$CC$subset$go,
                                     pathways), T, F)

openxlsx::write.xlsx(COMP.2DN_vs_Tumor,
                     file.path(date, results_dir, "2DN_vs_Tumor_GO_ACTIVATION.xlsx"))

## 2D hypoxia vs Tumor
COMP.2DH_vs_Tumor <- bind_df(SIMMAT.2DH_vs_Tumor, KEGG.2DH_vs_Tumor$df)
COMP.2DH_vs_Tumor$category <- factor(COMP.2DH_vs_Tumor$category, 
                                     levels=c("BP", "MF", "CC", "KEGG"))

COMP.2DH_vs_Tumor <- make_GObase(COMP.2DH_vs_Tumor, results.2DH_vs_Tumor$sig_df)
COMP.2DH_vs_Tumor$slim <- ifelse(COMP.2DH_vs_Tumor$id %in% 
                                   c(processes,
                                     REDUCED.2DH_vs_Tumor$MF$subset$go,
                                     REDUCED.2DH_vs_Tumor$CC$subset$go,
                                     pathways), T, F)

openxlsx::write.xlsx(COMP.2DH_vs_Tumor,
                     file.path(date, results_dir, "2DH_vs_Tumor_GO_ACTIVATION.xlsx"))

## 3D organoids vs Tumor
COMP.3D_vs_Tumor <- bind_df(SIMMAT.3D_vs_Tumor, KEGG.3D_vs_Tumor$df)
COMP.3D_vs_Tumor$category <- factor(COMP.3D_vs_Tumor$category, 
                                     levels=c("BP", "MF", "CC", "KEGG"))

COMP.3D_vs_Tumor <- make_GObase(COMP.3D_vs_Tumor, results.3D_vs_Tumor$sig_df)
COMP.3D_vs_Tumor$slim <- ifelse(COMP.3D_vs_Tumor$id %in% 
                                   c(processes,
                                     REDUCED.3D_vs_Tumor$MF$subset$go,
                                     REDUCED.3D_vs_Tumor$CC$subset$go,
                                     pathways), T, F)

openxlsx::write.xlsx(COMP.3D_vs_Tumor,
                     file.path(date, results_dir, "3D_vs_Tumor_GO_ACTIVATION.xlsx"))

# -------------------------- In-vitro results ----------------------------------
# 1.1 Prepare gene sets for semantic similarity analysis
SIMMAT.2DH_vs_2DN <- list()
SIMMAT.2DH_vs_2DN <- GO.2DH_vs_2DN$df %>%
  dplyr::group_split(ONTOLOGY)

SIMMAT.2DH_vs_2DN <- setNames(SIMMAT.2DH_vs_2DN, ontology)

SIMMAT.3D_vs_2DN <- list()
SIMMAT.3D_vs_2DN <- GO.3D_vs_2DN$df %>%
  dplyr::group_split(ONTOLOGY)

SIMMAT.3D_vs_2DN <- setNames(SIMMAT.3D_vs_2DN, ontology)

# 1.2 Run the semantic similarity analysis
## 2D hypoxia vs 2D normoxia
SIM.2DH_vs_2DN <- make_simMatrix(SIMMAT.2DH_vs_2DN) # make similarity matrix
SCORES.2DH_vs_2DN <- lapply(SIMMAT.2DH_vs_2DN, function(x) setNames(-log10(x$pvalue), x$ID))
for(var in names(SCORES.2DH_vs_2DN)){
  SCORES.2DH_vs_2DN[[var]] <- SCORES.2DH_vs_2DN[[var]][which(names(SCORES.2DH_vs_2DN[[var]]) %in% row.names(SIM.2DH_vs_2DN[[var]]))]
}  # calculate the scores for the significant terms
REDUCED.2DH_vs_2DN <- list()
for(var in names(SIM.2DH_vs_2DN)){
  REDUCED.2DH_vs_2DN[[var]] <- get_reducedTerms(simm = SIM.2DH_vs_2DN[[var]], 
                                         scores = SCORES.2DH_vs_2DN[[var]], 
                                         0.9)
} # reduce collinearity

## 3D organoids vs 2D normoxia
SIM.3D_vs_2DN <- make_simMatrix(SIMMAT.3D_vs_2DN) # make similarity matrix
SCORES.3D_vs_2DN <- lapply(SIMMAT.3D_vs_2DN, function(x) setNames(-log10(x$pvalue), x$ID))
for(var in names(SCORES.3D_vs_2DN)){
  SCORES.3D_vs_2DN[[var]] <- SCORES.3D_vs_2DN[[var]][which(names(SCORES.3D_vs_2DN[[var]]) %in% row.names(SIM.3D_vs_2DN[[var]]) )]
} # calculate the scores for the significant terms
REDUCED.3D_vs_2DN <- list()
for(var in names(SIM.3D_vs_2DN)){
  REDUCED.3D_vs_2DN[[var]] <- get_reducedTerms(simm = SIM.3D_vs_2DN[[var]], 
                                         scores = SCORES.3D_vs_2DN[[var]], 
                                         0.9)
} # reduce collinearity

## 1.3 Create compound data frame for visualization
## 2D hypoxia vs 2D normoxia
COMP.2DH_vs_2DN <- bind_df(SIMMAT.2DH_vs_2DN, KEGG.2DH_vs_2DN$df)
COMP.2DH_vs_2DN$category <- factor(COMP.2DH_vs_2DN$category, 
                                   levels=c("BP", "MF", "CC", "KEGG"))

COMP.2DH_vs_2DN <- make_GObase(COMP.2DH_vs_2DN, results.2DH_vs_2DN$sig_df)
COMP.2DH_vs_2DN$slim <- ifelse(COMP.2DH_vs_2DN$id %in% 
                                   c(processes,
                                     REDUCED.2DH_vs_2DN$MF$subset$go,
                                     REDUCED.2DH_vs_2DN$CC$subset$go,
                                     pathways), T, F)

openxlsx::write.xlsx(COMP.2DH_vs_2DN, 
           file.path(date, results_dir, "2DH_vs_2DN_GO_ACTIVATION.xlsx"))

## 3D organoids vs 2D normoxia
COMP.3D_vs_2DN <- bind_df(SIMMAT.3D_vs_2DN, KEGG.3D_vs_2DN$df)
COMP.3D_vs_2DN$category <- factor(COMP.3D_vs_2DN$category, 
                                   levels=c("BP", "MF", "CC", "KEGG"))

COMP.3D_vs_2DN <- make_GObase(COMP.3D_vs_2DN, results.3D_vs_2DN$sig_df)
COMP.3D_vs_2DN$slim <- ifelse(COMP.3D_vs_2DN$id %in% 
                                   c(processes,
                                     REDUCED.3D_vs_2DN$MF$subset$go,
                                     REDUCED.3D_vs_2DN$CC$subset$go,
                                     pathways), T, F)

openxlsx::write.xlsx(COMP.3D_vs_2DN,
           file.path(date, results_dir, "3D_vs_2DN_GO_ACTIVATION.xlsx"))

####################################################
# 10.) Gene sets in interception                   #
####################################################
# ---------------------- Total results -----------------------------------------
# Create Venn diagram of shared differentially expressed genes 
simple_venn(results.2DN_vs_Tumor$sig_df, 
            results.2DH_vs_Tumor$sig_df, 
            results.3D_vs_Tumor$sig_df,
            "venn_diagram")

venn <- make_vennbase(results.2DN_vs_Tumor$sig_df, 
                        results.2DH_vs_Tumor$sig_df, 
                        results.3D_vs_Tumor$sig_df,
                        c("2DN", "2DH", "3D")) 

# get the membership of the genes in the regions
upset <- make_upsetbase(venn$table, c("2DN", "2DH", "3D"))
# get the centrum of the regions
arranged <- arrange_venn(upset, c("2DN", "2DH", "3D"),
                               extract_regions = T) 

set.seed(123)
# randomize the position of the genes in the regions
xy = rbind(
  calculateCircle(x = 0.00, y = 0.2886, r = 0.1, # 3-way intersection
                  noiseFun = function(x) (x + rnorm(1,0,0.1)), # add noise
                  steps = 376,randomDist = T, randomFun = rnorm), # 579 genes
  calculateEllipse(x = -0.7, y = -0.1266, a = 0.2, b = .2, angle = -240, 
                   noiseFun = function(x) (x + rnorm(1,0,0.1)),
                   steps = 83, randomDist = T, randomFun = rnorm),
  calculateEllipse(x = 0.7, y = -0.1266, a = 0.2, b = .2, angle = -120, 
                   noiseFun = function(x) (x + rnorm(1,0,0.1)),
                   steps = 18, randomDist = T, randomFun = rnorm),
  calculateEllipse(x = 0, y = 1.0992, a = 0.1, b = .2, angle = 0, 
                   noiseFun = function(x) (x + rnorm(1,0,0.1)),
                   steps = 216, randomDist = T, randomFun = rnorm),
  calculateCircle(x = 0.00, y = -1.2124, r = .3, 
                  noiseFun = function(x) (x + rnorm(1,0,0.2)),
                  steps = 317, randomDist = T, randomFun = rnorm),
  calculateCircle(x = 1.30, y = 1.0392, r = .3,  
                  noiseFun = function(x) (x + rnorm(1,0,0.2)),
                  steps = 58,randomDist = T, randomFun = rnorm),
  calculateCircle(x = -1.30, y = 1.0392, r = .3, 
                  noiseFun = function(x) (x + rnorm(1,0,0.2)),
                  steps = 58,randomDist = T, randomFun = rnorm)
)

upset <- upset %>%
  dplyr::mutate(region = dplyr::case_when(
    `2DN` & `2DH` & `3D` ~ "2DH-2DN-3D",
    `2DN` & !`2DH` & `3D` ~ "2DN-3D",
    !`2DN` & `2DH` & `3D` ~ "2DH-3D",
    `2DN` & `2DH` & !`3D` ~ "2DH-2DN",
    !`2DN` & !`2DH` & `3D` ~ "3D",
    !`2DN` & `2DH` & !`3D` ~ "2DH",
    `2DN` & !`2DH` & !`3D` ~ "2DN"
  ), # relabel the regions
  region = as.factor(region),
  x = xy[,1],
  y = xy[,2],
  ) # add the x and y coordinates of the genes

upset <- upset %>% 
  # add their log2FoldChange values to the data frame
  dplyr::mutate(log2FoldChange = rowMeans(upset[,5:7], na.rm = T)) %>% 
  # calculate mean log2FoldChange values for the intersection genes
  dplyr::arrange(desc(abs(log2FoldChange))*-1) 
# sort the genes by their log2FoldChange values

####################################################
# 11.) WGCNA analysis                              #
####################################################
expr.matrix <- merge.rec(list("IDs"= data.frame(ensemblID = results.2DN_vs_Tumor$df$ensemblID),
                              "2DN" = results.2DN_vs_Tumor$sig_df[,c("ensemblID","log2FoldChange")],
                              "2DH" = results.2DH_vs_Tumor$sig_df[,c("ensemblID","log2FoldChange")],
                              "3D" = results.3D_vs_Tumor$sig_df[,c("ensemblID","log2FoldChange")]),
                         by = "ensemblID", all = T) %>%
  tibble::column_to_rownames("ensemblID") %>%
  setNames(c("Tumor vs 2DN","Tumor vs 2DH","Tumor vs 3D"))

# 1.1 Prepare the data
library(WGCNA)
# Here we are not specifying a model, only using the model intersect
wgcna <- list()
wgcna$vst <- t(assay(vst(deseq$dds[rownames(expr.matrix),])))

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

wgcna$soft <- WGCNA::pickSoftThreshold(
  wgcna$vst,
  powerVector = powers,
  verbose = 5
)

wgcna$soft_df <- data.frame(wgcna$soft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(wgcna$soft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(wgcna$soft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()

temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
wgcna$net <- WGCNA::blockwiseModules(wgcna$vst,
                          # == Adjacency Function ==
                          power = 8,               
                          networkType = "signed",
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 15000,
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

readr::write_rds(wgcna$net,
                 file = file.path(date, results_dir, "wgcna_results.RDS")
)
# Convert labels to colors for plotting
mergedColors <- WGCNA::labels2colors(wgcna$net$colors)
# Plot the dendrogram and the module colors underneath
module_df <- data.frame(
  gene_id = names(wgcna$net$colors),
  colors = WGCNA::labels2colors(wgcna$net$colors)
)

svg(file.path(date, plots_dir, "dendrogram.svg"),
    width = 18, height = 10)
WGCNA::plotDendroAndColors(
  wgcna$net$dendrograms[[1]],
  mergedColors[wgcna$net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05)
dev.off()

lognorm.matrix <- make_matrix(assay(deseq$rld[rownames(expr.matrix),]), 
                              module_df[wgcna$net$blockGenes[[1]], "gene_id"])



dend <- hang.dendrogram(as.dendrogram(wgcna$net$dendrograms[[1]]), hang = 0.03)

heatmap.annot <- HeatmapAnnotation(
  patient = as.factor(coldata$patient),
  condition = as.factor(coldata$condition),
  show_annotation_name = F,
  col =list(patient = c("593" = "navy", "673" = "khaki"),
            condition = c("normoxia.2D" = "steelblue",
                          "hypoxia.2D" = "salmon",
                          "physioxia.3D" = "red",
                          "physioxia.Tumour" = "darkred")))

modul.annot <- rowAnnotation(
  module = as.factor(module_df[wgcna$net$blockGenes[[1]], "colors"]),
  col = list(module = c("black"="black","blue"="blue","brown"="brown",
                        "green"="green","grey"="grey","magenta"="magenta",
                        "pink"="pink","red"="red","turquoise"="turquoise","yellow"="yellow")),
  show_annotation_name = T, annotation_label = "Module colors", 
  width = unit(1, "in"))
  
hm <- Heatmap(lognorm.matrix, name = "Normalized expression",
              cluster_rows = dend,
              show_column_dend = T, show_column_names = F,
              show_row_names = F, show_row_dend = F,
              column_dend_reorder = T,
              col = colorRamp2(c(min(lognorm.matrix), 0, max(lognorm.matrix)), 
                               c("green","black","red")),
              left_annotation = modul.annot,
              bottom_annotation = heatmap.annot)

exp <- Heatmap(as.matrix(expr.matrix), name = "log2FC",
               na_col = "white", border = T,
               column_split = 3,
               cluster_rows = dend,
               column_title = "Expression change",
               show_column_names = T, show_row_names = F)
  

svg(file.path(date, plots_dir,"total_heatmap.svg"),
    width = 18, height = 10)
draw(hm + exp, merge_legend = T, padding = unit(c(15, 2, 2, 2), "mm"))
dev.off()

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(wgcna$vst, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$condition = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  tidyr::pivot_longer(-condition) %>%
  dplyr::mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

module_plot <- mME %>% ggplot(., aes(x=condition, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

ggsave(file.path(date, plots_dir, "module_trait_relationships.png"),
       module_plot, width = 12, height = 10)

dt <- mME %>% 
  tidyr::pivot_wider(names_from = name, values_from = value) %>% 
  tibble::column_to_rownames("condition")

mod_ht <- draw(ComplexHeatmap::Heatmap(
  t(dt),
  name = "Module eigengenes",
  col = colorRamp2(c(min(dt), 0, max(dt)),
                   c("blue", "white", "red")),
  show_column_dend = T, show_column_names = T,
  show_row_names = T, show_row_dend = T,
  column_dend_reorder = T, row_split = 3,
  row_dend_reorder = T,
  clustering_distance_columns = "euclidean",
  column_title = "Module eigengenes",
  row_title = "Modules",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.3f", t(dt)[i, j]), x, y, gp = gpar(fontsize = 10))},
  row_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 12),
  row_title_gp = gpar(fontsize = 12),
  bottom_annotation = heatmap.annot,
  heatmap_legend_param = list(title = "Module eigengenes"),
  width = unit(10, "cm"),
  height = unit(10, "cm")
), merge_legend = TRUE)

# pick out a few modules of interest here
modules_of_interest = c("blue", "grey")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get the expression data for the genes in the module
submod_df = data.frame(t(wgcna$vst)) %>%
  dplyr::filter(row.names(.) %in% submod$gene_id) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  tidyr::pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors,
    module = as.factor(module)
  )


