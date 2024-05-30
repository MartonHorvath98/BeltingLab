# ---------------------------------------------------------------------------- #
# -            SET UP WORKING DIRECTORY AND DIRECTORY TREE                   - #
# ---------------------------------------------------------------------------- #
# 1.) Set downstream path
## Create the sub folders for: results, data, and pictures
data <- "data"
if (!dir.exists(file.path(data))) {
  dir.create(file.path(data)) # create the data folder
} 
plots_dir <- "plots" # plots directory
results_dir <- "results" # results directory
date <- format(Sys.Date(), "%Y-%m-%d") # get the current date
if (!dir.exists(file.path(date))) {
  dir.create(file.path(date)) # create the dated results folder
  dir.create(file.path(date, results_dir)) # create the results folder
  dir.create(file.path(date, plots_dir)) # create the plots folder
}

# 2.) Load functions for the analyses and required packages
source("functions.R")
source("packages.R")
source("helpers.R")

# 3.) Load analysis data              
if (exists("readcounts") == F) {
  # copy the read counts files generated with featureCounts to the data folder
  mrna.path <- file.path("../../results/06_counts/")
  mrna.files <- list.files(mrna.path, pattern = "counts.txt$", full.names = T)
  file.copy(from = mrna.files, to = file.path(data))
  
  # input is a tab-delimited csv, with quantification parameters in the first row 
  mrna.reads <- lapply(list.files(file.path(data), 
                                  pattern = "counts.txt$",
                                  full.names = T), 
                       read.csv, header = T, sep = "\t", skip = 1)
  
  # feature names are in the first column, counts in the 7th, the rest are gene length,
  # orientation, etc. We only need the gene ID and the counts:
  mrna.counts <- lapply(mrna.reads, 
                        function(x) { x  %>%
                            dplyr::select(c(1,7)) %>%
                            setNames(c("Geneid","Counts")) %>% 
                            dplyr::mutate(Geneid = stringr::str_remove(Geneid, "gene:"))}
  )
  # Recursively merge the list of count matrices into a single data frame
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

# 4.) Load count read counts matrix
readcounts <- read.csv(file = file.path(data,"readcounts.csv"),
                       sep = ",", header = T, na.strings = NA, row.names = 1)
readcounts <- as.matrix(readcounts)

# 5.) Prepare metadata table based on the sample names
samples <- colnames(readcounts)
coldata <- data.frame("samplenames" = samples) %>%
  # extract information from the sample names
  dplyr::mutate(samplenames = as.factor(samplenames)) %>% 
  dplyr::mutate(patient = dplyr::case_when( # extract the patient IDs
    stringr::str_detect(samples, "593") ~ "593",
    stringr::str_detect(samples, "673") ~ "673"
    ),
    patient = factor(patient, levels = c("593","673"))) %>%
  dplyr::mutate(dimension = dplyr::case_when( # extract the dimension
    stringr::str_detect(samples, ".2D") ~ "2D",
    stringr::str_detect(samples, ".3D") ~ "3D",
    TRUE ~ "Tumour"
    ),
    dimension = factor(dimension, levels = c("Tumour","2D","3D"))) %>% 
  dplyr::mutate(oxygen = dplyr::case_when( # extract the growth conditions
    stringr::str_detect(samples, "H.") ~ "hypoxia",
    stringr::str_detect(samples, "N.") ~ "normoxia",
    TRUE ~ "physioxia"
    ),
    oxygen = factor(oxygen, levels = c("physioxia","normoxia","hypoxia"))) %>%
  dplyr::mutate(run = dplyr::case_when( # extract the technical replicates
    stringr::str_detect(samples, ".1$") ~ "run1",
    stringr::str_detect(samples, ".2$") ~ "run2",
    stringr::str_detect(samples, ".3$") ~ "run3"
    ),
    run = factor(run, levels = c("run1","run2","run3"))) %>%
  dplyr::mutate(
    condition = factor( # set up the factor describing the experimental design
      str_glue("{oxygen}.{dimension}"),
      levels = c("normoxia.2D","hypoxia.2D","physioxia.3D","physioxia.Tumour")))

# ---------------------------------------------------------------------------- #
#   GENOME-WIDE ANALYSIS OF THE EFFECT OF OXYGEN AND DIMENSION ON GENE EXPR.   #
# ---------------------------------------------------------------------------- #

# ------- In vitro results: using 2D normoxia as baseline expression --------- #

# 1.) Differential gene expression analysis
## Set up the metadata table
coldata.invitro <- coldata
## Run DESeq2
deseq.invitro <- make_deseq(matrix = readcounts, # expression matrix
                            coldata = coldata.invitro, # metadata table
                            design = "condition + condition:patient") # design formula
resultsNames(deseq.invitro$dds)

## Extract significant results
results.2DH_vs_2DN <- get_results(dds = deseq.invitro$dds, # DESeq2 object
                                  sig_log2FC = 1.5, sig_pval = 0.05, # significance thresholds
                                  contrast=list("condition_hypoxia.2D_vs_normoxia.2D"),
                                  name = "condition_hypoxia.2D_vs_normoxia.2D")
# Number of up- and down-regulated genes between 2D normoxia vs 2D hypoxia
table(results.2DH_vs_2DN$sig_df$significance)
# Signif. down-regulated   Signif. up-regulated 
# 7                        123 

results.3D_vs_2DN <- get_results(dds = deseq.invitro$dds, # DESeq2 object
                                 sig_log2FC = 1.5, sig_pval = 0.05, # significance thresholds
                                 contrast=list("condition_physioxia.3D_vs_normoxia.2D"),
                                 name = "condition_physioxia.3D_vs_normoxia.2D")
# Number of up- and down-regulated genes between 3D organoids vs 2D normoxia
table(results.3D_vs_2DN$sig_df$significance)
# Signif. down-regulated   Signif. up-regulated 
# 234                      160

results.Tumor_vs_2DN <- get_results(dds = deseq.invitro$dds, # DESeq2 object
                                    sig_log2FC = 1.5, sig_pval = 0.05, # significance thresholds
                                    contrast=list("condition_physioxia.Tumour_vs_normoxia.2D"),
                                    name = "condition_physioxia.Tumour_vs_normoxia.2D")
# Number of up- and down-regulated genes between Tumor vs 2D normoxia
table(results.Tumor_vs_2DN$sig_df$significance)
# Signif. down-regulated   Signif. up-regulated 
# 296                      436

## Save results to excel files
sapply(c("df","sig_df"), function(x){
  # 2D hypoxia vs 2D normoxia
  write.csv(results.2DH_vs_2DN[[x]],
            file.path(date, results_dir, paste0("2DH_vs_2DN_",x,".csv")),
            row.names = T)
  # 3D organoids vs 2D normoxia
  write.csv(results.3D_vs_2DN[[x]], 
            file.path(date, results_dir, paste0("3D_vs_2DN_",x,".csv")),
            row.names = T)
  # Tumour vs 2D normoxia
  write.csv(results.Tumor_vs_2DN[[x]], 
            file.path(date, results_dir, paste0("Tumor_vs_2DN_",x,".csv")),
            row.names = T)
})

# 2.) Over-representation (ORA) and Gene Set Enrichment Analysis (GSEA)
## KEGG pathways
## extract entrez IDs for gene set of interest and background
genes.2DH_vs_2DN <- get_genelist(results.2DH_vs_2DN$df, # 2D hypoxia vs 2D normoxia
                                 filter = results.2DH_vs_2DN$sig_df$entrezID)
genes.3D_vs_2DN <- get_genelist(results.3D_vs_2DN$df, # 3D organoids vs 2D normoxia
                                filter = results.3D_vs_2DN$sig_df$entrezID)
genes.Tumor_vs_2DN <- get_genelist(results.Tumor_vs_2DN$df, # Tumour vs 2D normoxia
                                  filter = results.Tumor_vs_2DN$sig_df$entrezID)

## Run the ORA analysis
KEGG.2DH_vs_2DN <- list()
KEGG.2DH_vs_2DN$df <- kegg_results(genes.2DH_vs_2DN$interest, # gene set of interest
                                   genes.2DH_vs_2DN$background) # background gene set

KEGG.3D_vs_2DN <- list()
KEGG.3D_vs_2DN$df <- kegg_results(genes.3D_vs_2DN$interest, # gene set of interest
                                  genes.3D_vs_2DN$background) # background gene set

KEGG.Tumor_vs_2DN <- list()
KEGG.Tumor_vs_2DN$df <- kegg_results(genes.Tumor_vs_2DN$interest, # gene set of interest
                                     genes.Tumor_vs_2DN$background) # background gene set
## Save the results
openxlsx::write.xlsx(KEGG.2DH_vs_2DN$df, 
                     file.path(date, results_dir, "2DH_vs_2DN_KEGG.xlsx"))

openxlsx::write.xlsx(KEGG.3D_vs_2DN$df,
                     file.path(date, results_dir, "3D_vs_2DN_KEGG.xlsx"))

openxlsx::write.xlsx(KEGG.Tumor_vs_2DN$df,
                     file.path(date, results_dir, "Tumor_vs_2DN_KEGG.xlsx"))

## GO terms
## run the GSEA analysis
GO.2DH_vs_2DN <- list()
GO.2DH_vs_2DN$df <- go_results(genes.2DH_vs_2DN$interest, # gene set of interest
                               genes.2DH_vs_2DN$background, # background gene set
                               type = 'ALL')
GO.3D_vs_2DN <- list()
GO.3D_vs_2DN$df <- go_results(genes.3D_vs_2DN$interest, # gene set of interest
                              genes.2DH_vs_Tumor$background, # background gene set
                              type = 'ALL')

GO.Tumor_vs_2DN <- list()
GO.Tumor_vs_2DN$df <- go_results(genes.Tumor_vs_2DN$interest, # gene set of interest
                                genes.Tumor_vs_2DN$background, # background gene set
                                type = 'ALL')

## Save the results
openxlsx::write.xlsx(GO.2DH_vs_2DN$df,
                     file.path(date, results_dir, "2DH_vs_2DN_GO.xlsx"))

openxlsx::write.xlsx(GO.3D_vs_2DN$df,
                     file.path(date, results_dir, "3D_vs_2DN_GO.xlsx"))

openxlsx::write.xlsx(GO.Tumor_vs_2DN$df,
                     file.path(date, results_dir, "Tumor_vs_2DN_GO.xlsx"))

# 3.) Semantic similarity analysis                 
## Prepare gene sets for semantic similarity analysis
SIMMAT.2DH_vs_2DN <- list()
SIMMAT.2DH_vs_2DN <- GO.2DH_vs_2DN$df %>%
  dplyr::group_split(ONTOLOGY) # split data frame by category
SIMMAT.2DH_vs_2DN <- setNames(SIMMAT.2DH_vs_2DN, 
                              levels(GO.2DH_vs_2DN$df$ONTOLOGY))
SIMMAT.3D_vs_2DN <- list()
SIMMAT.3D_vs_2DN <- GO.3D_vs_2DN$df %>%
  dplyr::group_split(ONTOLOGY)
SIMMAT.3D_vs_2DN <- setNames(SIMMAT.3D_vs_2DN, 
                             levels(GO.3D_vs_2DN$df$ONTOLOGY))
SIMMAT.Tumor_vs_2DN <- list()
SIMMAT.Tumor_vs_2DN <- GO.Tumor_vs_2DN$df %>%
  dplyr::group_split(ONTOLOGY)
SIMMAT.Tumor_vs_2DN <- setNames(SIMMAT.Tumor_vs_2DN, 
                                levels(GO.Tumor_vs_2DN$df$ONTOLOGY))

## Run the semantic similarity analysis
SIM.2DH_vs_2DN <- make_simMatrix(SIMMAT.2DH_vs_2DN) # make similarity matrix
SCORES.2DH_vs_2DN <- lapply(SIMMAT.2DH_vs_2DN, # calculate the scores based on p-values
                            function(x) setNames(-log10(x$pvalue), x$ID))
for(var in names(SCORES.2DH_vs_2DN)){
  SCORES.2DH_vs_2DN[[var]] <- SCORES.2DH_vs_2DN[[var]][which(names(SCORES.2DH_vs_2DN[[var]]) %in% row.names(SIM.2DH_vs_2DN[[var]]))]
}  
SIM.3D_vs_2DN <- make_simMatrix(SIMMAT.3D_vs_2DN) # make similarity matrix
SCORES.3D_vs_2DN <- lapply(SIMMAT.3D_vs_2DN, # calculate the scores based on p-values
                           function(x) setNames(-log10(x$pvalue), x$ID))
for(var in names(SCORES.3D_vs_2DN)){
  SCORES.3D_vs_2DN[[var]] <- SCORES.3D_vs_2DN[[var]][which(names(SCORES.3D_vs_2DN[[var]]) %in% row.names(SIM.3D_vs_2DN[[var]]) )]
}
SIM.Tumor_vs_2DN <- make_simMatrix(SIMMAT.Tumor_vs_2DN) # make similarity matrix
SCORES.Tumor_vs_2DN <- lapply(SIMMAT.Tumor_vs_2DN, # calculate the scores based on p-values
                              function(x) setNames(-log10(x$pvalue), x$ID))
for(var in names(SCORES.Tumor_vs_2DN)){
  SCORES.Tumor_vs_2DN[[var]] <- SCORES.Tumor_vs_2DN[[var]][which(names(SCORES.Tumor_vs_2DN[[var]]) %in% row.names(SIM.Tumor_vs_2DN[[var]]) )]
}

## Reduce collinearity
REDUCED.2DH_vs_2DN <- list()
for(var in names(SIM.2DH_vs_2DN)){
  # similarity threshold is set to 'high' = 0.9
  REDUCED.2DH_vs_2DN[[var]] <- get_reducedTerms(simm = SIM.2DH_vs_2DN[[var]], 
                                                scores = SCORES.2DH_vs_2DN[[var]], 
                                                limit = 1) # select top hit per cluster
}
REDUCED.3D_vs_2DN <- list()
for(var in names(SIM.3D_vs_2DN)){
  # similarity threshold is set to 'high' = 0.9
  REDUCED.3D_vs_2DN[[var]] <- get_reducedTerms(simm = SIM.3D_vs_2DN[[var]], 
                                               scores = SCORES.3D_vs_2DN[[var]], 
                                               limit = 1) # select top hit per cluster
}
REDUCED.Tumor_vs_2DN <- list()
for(var in names(SIM.Tumor_vs_2DN)){
  # similarity threshold is set to 'high' = 0.9
  REDUCED.Tumor_vs_2DN[[var]] <- get_reducedTerms(simm = SIM.Tumor_vs_2DN[[var]], 
                                                  scores = SCORES.Tumor_vs_2DN[[var]], 
                                                  limit = 1) # select top hit per cluster
}

## Merge KEGG and GO results and add activation score
COMP.2DH_vs_2DN <- bind_df(SIMMAT.2DH_vs_2DN, KEGG.2DH_vs_2DN$df)
COMP.2DH_vs_2DN <- make_GObase(COMP.2DH_vs_2DN, results.2DH_vs_2DN$sig_df)
COMP.2DH_vs_2DN <- COMP.2DH_vs_2DN %>% 
  dplyr::mutate(slim = case_when(
    ID %in% c(REDUCED.2DH_vs_2DN$BP$subset$go,
              REDUCED.2DH_vs_2DN$MF$subset$go,
              REDUCED.2DH_vs_2DN$CC$subset$go,
              pathways) & 
    abs(zscore) >= 2 ~ TRUE,
    TRUE ~ FALSE))
COMP.3D_vs_2DN <- bind_df(SIMMAT.3D_vs_2DN, KEGG.3D_vs_2DN$df)
COMP.3D_vs_2DN <- make_GObase(COMP.3D_vs_2DN, results.3D_vs_2DN$sig_df)
COMP.3D_vs_2DN <- COMP.3D_vs_2DN %>% 
  dplyr::mutate(slim = case_when(
    ID %in% c(REDUCED.3D_vs_2DN$BP$subset$go,
              REDUCED.3D_vs_2DN$MF$subset$go,
              REDUCED.3D_vs_2DN$CC$subset$go,
              pathways) & 
    abs(zscore) >= 2 ~ TRUE,
    TRUE ~ FALSE))
COMP.Tumor_vs_2DN <- bind_df(SIMMAT.Tumor_vs_2DN, KEGG.Tumor_vs_2DN$df)
COMP.Tumor_vs_2DN <- make_GObase(COMP.Tumor_vs_2DN, results.Tumor_vs_2DN$sig_df)
COMP.Tumor_vs_2DN <- COMP.Tumor_vs_2DN %>% 
  dplyr::mutate(slim = case_when(
    ID %in% c(REDUCED.Tumor_vs_2DN$BP$subset$go,
              REDUCED.Tumor_vs_2DN$MF$subset$go,
              REDUCED.Tumor_vs_2DN$CC$subset$go,
              pathways) & 
    abs(zscore) >= 2 ~ TRUE,
    TRUE ~ FALSE))

## Save the results
openxlsx::write.xlsx(COMP.2DH_vs_2DN,
                     file.path(date, results_dir, "2DH_vs_2DN_ACTIVATION.xlsx"))
openxlsx::write.xlsx(COMP.3D_vs_2DN,
                     file.path(date, results_dir, "3D_vs_2DN_ACTIVATION.xlsx"))
openxlsx::write.xlsx(COMP.Tumor_vs_2DN,
                     file.path(date, results_dir, "Tumor_vs_2DN_ACTIVATION.xlsx"))

# -------- Total results: using tumor samples as baseline expression --------- #

# 1.) Differential gene expression analysis
## Relevel the coldata to have tumor expression as reference
coldata$condition <- relevel(coldata$condition, ref = "physioxia.Tumour")
## Run DESeq2
deseq <- make_deseq(matrix = readcounts, # expression matrix
                    coldata = coldata, # metadata table
                    design = "condition + condition:patient") # design formula
resultsNames(deseq$dds)

## Extract significant results
results.2DN_vs_Tumor <- get_results(dds = deseq$dds, # DESeq2 object
                                    sig_log2FC = 1.5, sig_pval = 0.05, # significance thresholds
                                    contrast=list("condition_normoxia.2D_vs_physioxia.Tumour"),
                                    name = "condition_normoxia.2D_vs_physioxia.Tumour")
# Number of up- and down-regulated genes between 2D normoxia and Tumor samples
table(results.2DN_vs_Tumor$sig_df$significance)
# Signif. down-regulated   Signif. up-regulated 
# 436                      297 

results.2DH_vs_Tumor <- get_results(dds = deseq$dds, # DESeq2 object
                                    sig_log2FC = 1.5, sig_pval = 0.05, # significance thresholds
                                    contrast=list("condition_hypoxia.2D_vs_physioxia.Tumour"),
                                    name = "condition_hypoxia.2D_vs_physioxia.Tumour")
# Number of up- and down-regulated genes between 2D hypoxia and Tumor samples
table(results.2DH_vs_Tumor$sig_df$significance)
# Signif. down-regulated   Signif. up-regulated 
# 386                      300 

results.3D_vs_Tumor <- get_results(dds = deseq$dds, # DESeq2 object
                                   sig_log2FC = 1.5, sig_pval = 0.05, # significance thresholds
                                   contrast=list("condition_physioxia.3D_vs_physioxia.Tumour"),
                                   name = "condition_physioxia.3D_vs_physioxia.Tumour")
# Number of up- and down-regulated genes between 3D organoids and Tumor samples
table(results.3D_vs_Tumor$sig_df$significance)
# Signif. down-regulated   Signif. up-regulated
# 508                      286

## Save results to excel files
sapply(c("df","sig_df"), function(x){
  # Tumor samples vs 2D normoxia
  write.csv(results.2DN_vs_Tumor[[x]],
            file.path(date, results_dir, paste0("2DN_vs_Tumor_",x,".csv")),
            row.names = T)
  # Tumor samples vs 2D hypoxia
  write.csv(results.2DH_vs_Tumor[[x]], 
            file.path(date, results_dir, paste0("2DH_vs_Tumor_",x,".csv")),
            row.names = T)
  # Tumour samples vs 3D organoids
  write.csv(results.3D_vs_Tumor[[x]], 
            file.path(date, results_dir, paste0("3D_vs_Tumor_",x,".csv")),
            row.names = T)
})

# 2.) Over-representation (ORA) and Gene Set Enrichment Analysis (GSEA)
# KEGG pathways 
## extract entrez IDs for the gene set of interest and background
genes.2DN_vs_Tumor <- get_genelist(results.2DN_vs_Tumor$df, # 2D normoxia vs Tumour
                                   filter = results.2DN_vs_Tumor$sig_df$entrezID)
genes.2DH_vs_Tumor <- get_genelist(results.2DH_vs_Tumor$df, # 2D hypoxia vs Tumour
                                   filter = results.2DH_vs_Tumor$sig_df$entrezID)
genes.3D_vs_Tumor <- get_genelist(results.3D_vs_Tumor$df, # 3D organoids vs Tumour
                                  filter = results.3D_vs_Tumor$sig_df$entrezID)

## Run the ORA analysis
KEGG.2DN_vs_Tumor <- list()
KEGG.2DN_vs_Tumor$df <- kegg_results(genes.2DN_vs_Tumor$interest,
                                     genes.2DN_vs_Tumor$background)

KEGG.2DH_vs_Tumor <- list()
KEGG.2DH_vs_Tumor$df <- kegg_results(genes.2DH_vs_Tumor$interest,
                                     genes.2DH_vs_Tumor$background)

KEGG.3D_vs_Tumor <- list()
KEGG.3D_vs_Tumor$df <- kegg_results(genes.3D_vs_Tumor$interest,
                                    genes.3D_vs_Tumor$background)

## Save the results
openxlsx::write.xlsx(KEGG.2DN_vs_Tumor$df, 
                     file.path(date, results_dir, "2DN_vs_Tumor_KEGG.xlsx"))

openxlsx::write.xlsx(KEGG.2DH_vs_Tumor$df,
                     file.path(date, results_dir, "2DH_vs_Tumor_KEGG.xlsx"))

openxlsx::write.xlsx(KEGG.3D_vs_Tumor$df,
                     file.path(date, results_dir, "3D_vs_Tumor_KEGG.xlsx"))

## GO terms
## run the GSEA analysis
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

## Save the results
openxlsx::write.xlsx(GO.2DN_vs_Tumor$df,
                     file.path(date, results_dir, "2DN_vs_Tumor_GO.xlsx"))

openxlsx::write.xlsx(GO.2DH_vs_Tumor$df,
                     file.path(date, results_dir, "2DH_vs_Tumor_GO.xlsx"))

openxlsx::write.xlsx(GO.3D_vs_Tumor$df,
                     file.path(date, results_dir, "3D_vs_Tumor_GO.xlsx"))

# 3.) Semantic similarity analysis
## Prepare gene sets for semantic similarity analysis
SIMMAT.2DN_vs_Tumor <- list()
SIMMAT.2DN_vs_Tumor <- GO.2DN_vs_Tumor$df %>%
  dplyr::group_split(ONTOLOGY)
SIMMAT.2DN_vs_Tumor <- setNames(SIMMAT.2DN_vs_Tumor, 
                                levels(GO.2DN_vs_Tumor$df$ONTOLOGY))

SIMMAT.2DH_vs_Tumor <- list()
SIMMAT.2DH_vs_Tumor <- GO.2DH_vs_Tumor$df %>%
  dplyr::group_split(ONTOLOGY)
SIMMAT.2DH_vs_Tumor <- setNames(SIMMAT.2DH_vs_Tumor, 
                                levels(GO.2DH_vs_Tumor$df$ONTOLOGY))

SIMMAT.3D_vs_Tumor <- list()
SIMMAT.3D_vs_Tumor <- GO.3D_vs_Tumor$df %>%
  dplyr::group_split(ONTOLOGY)
SIMMAT.3D_vs_Tumor <- setNames(SIMMAT.3D_vs_Tumor, 
                               levels(GO.3D_vs_Tumor$df$ONTOLOGY))

## Run the semantic similarity analysis
SIM.2DN_vs_Tumor <- make_simMatrix(SIMMAT.2DN_vs_Tumor) # make similarity matrix
SCORES.2DN_vs_Tumor <- lapply(SIMMAT.2DN_vs_Tumor, # calculate the scores based on p-values
                              function(x) setNames(-log10(x$pvalue), x$ID))
for(var in names(SCORES.2DN_vs_Tumor)){
  SCORES.2DN_vs_Tumor[[var]] <- SCORES.2DN_vs_Tumor[[var]][which(names(SCORES.2DN_vs_Tumor[[var]]) %in% row.names(SIM.2DN_vs_Tumor[[var]]))]
}  
SIM.2DH_vs_Tumor <- make_simMatrix(SIMMAT.2DH_vs_Tumor) # make similarity matrix
SCORES.2DH_vs_Tumor <- lapply(SIMMAT.2DH_vs_Tumor, # calculate the scores based on p-values
                              function(x) setNames(-log10(x$pvalue), x$ID))
for (var in names(SCORES.2DH_vs_Tumor)){
  SCORES.2DH_vs_Tumor[[var]] <- SCORES.2DH_vs_Tumor[[var]][which(names(SCORES.2DH_vs_Tumor[[var]]) %in% row.names(SIM.2DH_vs_Tumor[[var]]))]
} 
SIM.3D_vs_Tumor <- make_simMatrix(SIMMAT.3D_vs_Tumor) # make similarity matrix
SCORES.3D_vs_Tumor <- lapply(SIMMAT.3D_vs_Tumor, # calculate the scores based on p-values
                             function(x) setNames(-log10(x$pvalue), x$ID))
for(var in names(SCORES.3D_vs_Tumor)){
  SCORES.3D_vs_Tumor[[var]] <- SCORES.3D_vs_Tumor[[var]][which(names(SCORES.3D_vs_Tumor[[var]]) %in% row.names(SIM.3D_vs_Tumor[[var]]) )]
} 

## Reduce collinearity
REDUCED.2DN_vs_Tumor <- list()
for(var in names(SIM.2DH_vs_Tumor)){
  # similarity threshold is set to 'high' = 0.9
  REDUCED.2DN_vs_Tumor[[var]] <- get_reducedTerms(simm = SIM.2DN_vs_Tumor[[var]], 
                                                  scores = SCORES.2DN_vs_Tumor[[var]], 
                                                  limit = 1) # select top hit per cluster
}
REDUCED.2DH_vs_Tumor <- list()
for(var in names(SIM.2DH_vs_Tumor)){
  # similarity threshold is set to 'high' = 0.9
  REDUCED.2DH_vs_Tumor[[var]] <- get_reducedTerms(simm = SIM.2DH_vs_Tumor[[var]], 
                                                  scores = SCORES.2DH_vs_Tumor[[var]],
                                                  limit = 1) # select top hit per cluster
}
REDUCED.3D_vs_Tumor <- list()
for(var in names(SIM.3D_vs_Tumor)){
  REDUCED.3D_vs_Tumor[[var]] <- get_reducedTerms(simm = SIM.3D_vs_Tumor[[var]], 
                                                 scores = SCORES.3D_vs_Tumor[[var]], 
                                                 limit = 1) # select top hit per cluster
}

## Merge KEGG and GO results and add activation score
COMP.2DN_vs_Tumor <- bind_df(SIMMAT.2DN_vs_Tumor, KEGG.2DN_vs_Tumor$df)
COMP.2DN_vs_Tumor <- make_GObase(COMP.2DN_vs_Tumor, results.2DN_vs_Tumor$sig_df)
COMP.2DN_vs_Tumor <- COMP.Tumor_vs_2DN %>% 
  dplyr::mutate(slim = case_when(
    ID %in% c(REDUCED.2DN_vs_Tumor$BP$subset$go,
              REDUCED.2DN_vs_Tumor$MF$subset$go,
              REDUCED.2DN_vs_Tumor$CC$subset$go,
              pathways) & 
      abs(zscore) >= 2 ~ TRUE,
    TRUE ~ FALSE))
COMP.2DH_vs_Tumor <- bind_df(SIMMAT.2DH_vs_Tumor, KEGG.2DH_vs_Tumor$df)
COMP.2DH_vs_Tumor <- make_GObase(COMP.2DH_vs_Tumor, results.2DH_vs_Tumor$sig_df)
COMP.2DH_vs_Tumor <- COMP.2DH_vs_Tumor %>% 
  dplyr::mutate(slim = case_when(
    ID %in% c(REDUCED.2DH_vs_Tumor$BP$subset$go,
              REDUCED.2DH_vs_Tumor$MF$subset$go,
              REDUCED.2DH_vs_Tumor$CC$subset$go,
              pathways) & 
      abs(zscore) >= 2 ~ TRUE,
    TRUE ~ FALSE))
COMP.3D_vs_Tumor <- bind_df(SIMMAT.3D_vs_Tumor, KEGG.3D_vs_Tumor$df)
COMP.3D_vs_Tumor <- make_GObase(COMP.3D_vs_Tumor, results.3D_vs_Tumor$sig_df)
COMP.3D_vs_Tumor <- COMP.3D_vs_Tumor %>% 
  dplyr::mutate(slim = case_when(
    ID %in% c(REDUCED.3D_vs_Tumor$BP$subset$go,
              REDUCED.3D_vs_Tumor$MF$subset$go,
              REDUCED.3D_vs_Tumor$CC$subset$go,
              pathways) & 
      abs(zscore) >= 2 ~ TRUE,
    TRUE ~ FALSE))

## Save the results
openxlsx::write.xlsx(COMP.2DN_vs_Tumor,
                     file.path(date, results_dir, "2DN_vs_Tumor_ACTIVATION.xlsx"))
openxlsx::write.xlsx(COMP.2DH_vs_Tumor,
                     file.path(date, results_dir, "2DH_vs_Tumor_ACTIVATION.xlsx"))
openxlsx::write.xlsx(COMP.3D_vs_Tumor,
                     file.path(date, results_dir, "3D_vs_Tumor_ACTIVATION.xlsx"))

# ---------------------------------------------------------------------------- #
# -         CLUSTER ANALYSIS OF CO-EXPRESSED GENES USING WGCNA               - #
# ---------------------------------------------------------------------------- #

# ----- Total results: using tumor samples as baseline expression ------------ #
# 1.) Prepare the data
expr.matrix <- merge.rec(
  list("IDs"= data.frame(ensemblID = results.2DN_vs_Tumor$df$ensemblID), # list of gene IDs
       "2DN" = results.2DN_vs_Tumor$sig_df[,c("ensemblID","log2FoldChange")],
       "2DH" = results.2DH_vs_Tumor$sig_df[,c("ensemblID","log2FoldChange")],
       "3D" = results.3D_vs_Tumor$sig_df[,c("ensemblID","log2FoldChange")]),
  by = "ensemblID", all = T) %>% # merge the data frames
  tibble::column_to_rownames("ensemblID") %>%
  setNames(c("2DN vs Tumor","2DH vs Tumor","3D vs Tumor"))

# 2.) Run WGCNA analysis
## Prepare the data: here we are not specifying a model, only using the model intersect
wgcna <- list()
wgcna <- make_wgcna(deseq = deseq$dds,
                    genes = expr.matrix,
                    initial = T) # Initial run to figure out soft power threshold
## Run the WGCNA analysis
wgcna$net <-  make_wgcna(deseq = deseq$dds,
                         genes =  expr.matrix, 
                         power = 8, 
                         initial = F)
readr::write_rds(wgcna$net,
                 file = file.path(date, results_dir, "wgcna_results.RDS")
)

## Module-trait relationships
wgcna$modules <- get_modules(wgcna$net, # get module genes end export dendogram
                             "WGCNA_colours_dendrogram")

wgcna$eigen <- moduleEigengenes(expr = wgcna$vst, # extract module eigen genes
                                colors = wgcna$modules$colors)$eigengenes

wgcna$connectivity <- signedKME(datExpr = wgcna$vst, # calculate connectivity 
                                datME = orderMEs(wgcna$eigen)) # within modules

wgcna$hubgenes <- apply(wgcna$connectivity, 2, function(x) { # export hub genes
  names(sort(x, decreasing = T)[1:10]) # top 10 genes with highest connectivity
  }) %>%                               # in each module
  as.data.frame() %>% 
  tidyr::pivot_longer(everything(), names_to = "module", values_to = "gene_id") %>% 
  dplyr::mutate(
    module = gsub("kME","",module),
    module = as.factor(module)
  )



# Module eigengene correlation with condition
MEs <- orderMEs(wgcna$eigen)
MEs <- MEs %>% 
  setNames(gsub("ME","",names(MEs)))

# Model matrix
MM <- cbind(table(coldata$samplenames, coldata$condition),
            table(coldata$samplenames, coldata$patient))

ME_corr <- WGCNA::cor(MEs, MM, use = "p")
ME_pval <- corPvalueStudent(ME_corr, ncol(MM))



AE <- make_matrix(t(wgcna$vst), wgcna$modules$gene_id) %>%
  as.data.frame() %>% 
  tibble::rownames_to_column("gene_id") %>%
  merge(.,  wgcna$modules, by="gene_id") %>% 
  dplyr::mutate(cluster = row_split[colors]) %>%
  dplyr::filter(! cluster %in% c("patient", "stable")) %>%
  tidyr::pivot_longer(-c("colors", "cluster", "gene_id"),
                      names_to = "samplenames", values_to = "averageExpr") %>% 
  merge(., coldata[,c("samplenames","condition")]) %>%
  dplyr::group_by(gene_id, cluster, colors, condition) %>% 
  dplyr::summarize(averageExpr = mean(averageExpr)) %>% 
  droplevels()

# 3.) ORA and GSE analyses on the clusters
## create a data frame with the module and cluster information
module_cluster <- module_df %>% 
  dplyr::filter(colors %in% c("black","yellow","red")) %>% 
  merge(., distinct(AE[,c("colors","cluster")]), by = "colors") %>% 
  dplyr::mutate(cluster = as.factor(cluster))

module_genes <- wgcna$hubgenes %>% 
  dplyr::filter(module %in% c("black","red","yellow")) %>% 
  dplyr::group_split(module) %>% 
  lapply(., function(x) {
    x %>% 
      dplyr::mutate(geneID =  AnnotationDbi::mapIds(org.Hs.eg.db, gene_id, 
                                                        keytype = "ENSEMBL", "SYMBOL")) %>% 
      dplyr::pull(geneID)
  }) %>% setNames(c("black","red","yellow"))

# 1.) Cluster1: 3D low genes ---------------------------------------------------
## n = 288
genes.3D <- get_genelist(results.3D_vs_2DN$df,
                         filter = module_cluster %>% 
                           dplyr::filter(cluster == "3D low") %>% 
                           dplyr::pull(gene_id) %>% 
                           AnnotationDbi::mapIds(org.Hs.eg.db, ., "ENTREZID", "ENSEMBL"))
# KEGG and GO analysis
KEGG.3D <- kegg_results(genes.3D$interest,
                        genes.3D$background)
write.xlsx(KEGG.3D, file.path(date, results_dir, "KEGG_3D.xlsx"))

GO.3D <- go_results(genes.3D$interest,
                    genes.3D$background,
                    type = 'ALL')
write.xlsx(GO.3D, file.path(date, results_dir, "GO_3D.xlsx"))

# 1.1 Prepare gene sets for semantic similarity analysis

SIMMAT.3D <- list()
SIMMAT.3D <- GO.3D %>%
  dplyr::group_split(ONTOLOGY)
SIMMAT.3D <- setNames(SIMMAT.3D, 
                      levels(GO.3D$ONTOLOGY))

# 1.2 Run the semantic similarity analysis

SIM.3D <- make_simMatrix(SIMMAT.3D) # make similarity matrix
SCORES.3D <- lapply(SIMMAT.3D, function(x) setNames(-log10(x$pvalue), x$ID))
for(var in names(SCORES.3D)){
  SCORES.3D[[var]] <- SCORES.3D[[var]][which(names(SCORES.3D[[var]]) %in% row.names(SIM.3D[[var]]))]
}  # calculate the scores for the significant terms
REDUCED.3D <- list()
for(var in names(SIM.3D)){
  REDUCED.3D[[var]] <- get_reducedTerms(simm = SIM.3D[[var]],
                                        scores = SCORES.3D[[var]], 
                                        limit = 1) # select top hit per cluster
} # reduce collinearity

## Merge KEGG and GO results and add activation score
COMP.3D <- bind_df(SIMMAT.3D, KEGG.3D)
COMP.3D <- make_GObase(COMP.3D, results.3D_vs_2DN$sig_df)
COMP.3D <- COMP.3D %>% 
  dplyr::mutate(slim = case_when(
    ID %in% c(REDUCED.3D$BP$subset$go,
              REDUCED.3D$MF$subset$go,
              REDUCED.3D$CC$subset$go,
              KEGG.3D$ID) &
      abs(zscore) >= 2 ~ TRUE,
    TRUE ~ FALSE))
openxlsx::write.xlsx(COMP.3D,
                     file.path(date, results_dir, "Module_3D_low_ACTIVATION.xlsx"))

# 2.) Cluster2: Tumor high genes -----------------------------------------------
## n = 498
genes.Tumor_high <- get_genelist(results.Tumor_vs_2DN$df,
                                 filter = module_cluster %>% 
                                   dplyr::filter(cluster == "Tumor high") %>% 
                                   dplyr::pull(gene_id) %>% 
                                   AnnotationDbi::mapIds(org.Hs.eg.db, ., "ENTREZID", "ENSEMBL"))
# KEGG and GO analysis
KEGG.Tumor_high <- kegg_results(genes.Tumor_high$interest,
                                genes.Tumor_high$background)
write.xlsx(KEGG.Tumor_high, file.path(date, results_dir, "KEGG_Tumor_high.xlsx"))

GO.Tumor_high <- go_results(genes.Tumor_high$interest,
                            genes.Tumor_high$background,
                            type = 'ALL')
write.xlsx(GO.Tumor_high, file.path(date, results_dir, "GO_Tumor_high.xlsx"))

# 1.1 Prepare gene sets for semantic similarity analysis

SIMMAT.Tumor_high <- list()
SIMMAT.Tumor_high <- GO.Tumor_high %>%
  dplyr::group_split(ONTOLOGY)
SIMMAT.Tumor_high <- setNames(SIMMAT.Tumor_high, 
                              levels(GO.Tumor_high$ONTOLOGY))

# 1.2 Run the semantic similarity analysis

SIM.Tumor_high <- make_simMatrix(SIMMAT.Tumor_high) # make similarity matrix
SCORES.Tumor_high <- lapply(SIMMAT.Tumor_high, function(x) setNames(-log10(x$pvalue), x$ID))
for(var in names(SCORES.Tumor_high)){
  SCORES.Tumor_high[[var]] <- SCORES.Tumor_high[[var]][which(names(SCORES.Tumor_high[[var]]) %in% row.names(SIM.Tumor_high[[var]]))]
}  # calculate the scores for the significant terms
REDUCED.Tumor_high <- list()
for(var in names(SIM.Tumor_high)){
  REDUCED.Tumor_high[[var]] <- get_reducedTerms(simm = SIM.Tumor_high[[var]],
                                                scores = SCORES.Tumor_high[[var]], 
                                                1)
} # reduce collinearity

## 1.3 Create compound data frame for visualization

COMP.Tumor_high <- bind_df(SIMMAT.Tumor_high, # combine KEGG and GO results
                           KEGG.Tumor_high)
COMP.Tumor_high <- make_GObase(COMP.Tumor_high, # calculate z-scores
                               results.Tumor_vs_2DN$sig_df)
COMP.Tumor_high <- COMP.Tumor_high %>% 
  dplyr::mutate(slim = case_when(
    ID %in% c(REDUCED.Tumor_high$BP$subset$go,
              REDUCED.Tumor_high$MF$subset$go,
              REDUCED.Tumor_high$CC$subset$go,
              KEGG.Tumor_high$ID) & 
    abs(zscore) >= 2 ~ TRUE,
    TRUE ~ FALSE))
openxlsx::write.xlsx(COMP.Tumor_high,
                     file.path(date, results_dir, "Module_Tumor_high_ACTIVATION.xlsx"))
# 3.) Cluster3: Tumor low genes ------------------------------------------------
## n = 1401
genes.Tumor_low <- get_genelist(results.Tumor_vs_2DN$df,
                                filter = module_cluster %>% 
                                  dplyr::filter(cluster == "Tumor low") %>% 
                                  dplyr::pull(gene_id) %>% 
                                  AnnotationDbi::mapIds(org.Hs.eg.db, ., "ENTREZID", "ENSEMBL"))

# KEGG and GO analysis
KEGG.Tumor_low <- kegg_results(genes.Tumor_low$interest,
                               genes.Tumor_low$background)
write.xlsx(KEGG.Tumor_low, file.path(date, results_dir, "KEGG_Tumor_low.xlsx"))

GO.Tumor_low <- go_results(genes.Tumor_low$interest,
                           genes.Tumor_low$background,
                           type = 'ALL')
write.xlsx(GO.Tumor_low, file.path(date, results_dir, "GO_Tumor_low.xlsx"))

# 1.1 Prepare gene sets for semantic similarity analysis

SIMMAT.Tumor_low <- list()
SIMMAT.Tumor_low <- GO.Tumor_low %>%
  dplyr::group_split(ONTOLOGY)
SIMMAT.Tumor_low <- setNames(SIMMAT.Tumor_low, 
                             levels(GO.Tumor_low$ONTOLOGY))

# 1.2 Run the semantic similarity analysis
SIM.Tumor_low <- make_simMatrix(SIMMAT.Tumor_low) # make similarity matrix
SCORES.Tumor_low <- lapply(SIMMAT.Tumor_low, function(x) setNames(-log10(x$pvalue), x$ID))
for(var in names(SCORES.Tumor_low)){
  SCORES.Tumor_low[[var]] <- SCORES.Tumor_low[[var]][which(names(SCORES.Tumor_low[[var]]) %in% row.names(SIM.Tumor_low[[var]]))]
}  # calculate the scores for the significant terms
REDUCED.Tumor_low <- list()
for(var in names(SIM.Tumor_low)){
  REDUCED.Tumor_low[[var]] <- get_reducedTerms(simm = SIM.Tumor_low[[var]],
                                               scores = SCORES.Tumor_low[[var]],
                                               1)
} # reduce collinearity

## 1.3 Create compound data frame for visualization
COMP.Tumor_low <- bind_df(SIMMAT.Tumor_low, # combine KEGG and GO results
                           KEGG.Tumor_low)
COMP.Tumor_low <- make_GObase(COMP.Tumor_low, # calculate z-scores
                              results.Tumor_vs_2DN$sig_df)
COMP.Tumor_low <- COMP.Tumor_low %>% 
  dplyr::mutate(slim = case_when(
    ID %in% c(REDUCED.Tumor_low$BP$subset$go,
              REDUCED.Tumor_low$CC$subset$go,
              pathways) & 
    abs(zscore) >= 2 ~ TRUE,
    TRUE ~ FALSE))
openxlsx::write.xlsx(COMP.Tumor_low,
                     file.path(date, results_dir, "Module_Tumor_low_ACTIVATION.xlsx"))

# 4.) Full, module activation   
## Combine the results
modules <- factor(c("3D low", "Tumor high", "Tumor low"))
  
COMP.full_modul <- list(
  "3D high" = COMP.3D,
  "Tumor_high" = COMP.Tumor_high,
  "Tumor_low" = COMP.Tumor_low)

COMP.full_modul <- lapply(modules, function(x){
  COMP.full_modul[[x]]$module <- x
  return(COMP.full_modul[[x]])
}) %>% 
  bind_rows()

# 5.) Deconvolute stromal cells                   
## Load the stromal cell type signatures
cell_type_rna <- list(
  "astrocytes" = read.table(file.path(data, "cell_astrocytes.tsv"),
                         header = T, sep = "\t"),
  "microglia" = read.table(file.path(data, "cell_microglia.tsv"),
                         header = T, sep = "\t"),
  "macrophages" = read.table(file.path(data, "cell_macrophage.tsv"),
                         header = T, sep = "\t"),
  "dendritic" = read.table(file.path(data, "cell_dendritic.tsv"),
                         header = T, sep = "\t"),
  "monocytes" = read.table(file.path(data, "cell_monocyte.tsv"),
                         header = T, sep = "\t"),
  "T-cells" = read.table(file.path(data, "cell_T-cell.tsv"),
                         header = T, sep = "\t"),
  "B-cells" = read.table(file.path(data, "cell_B-cell.tsv"),
                         header = T, sep = "\t")
)
cell_types <- names(cell_type_rna)
names(cell_types) <- c("Astrocytes", "Microglial cells",
                       "Macrophages", "dendritic cells",
                       "monocytes", "T-cells", "B-cells")

cell_type_rna <- lapply(cell_types, function(x){
  cells <- names(x)
  return(cell_type_rna[[x]] %>% 
           dplyr::mutate(
             TPM = str_extract(
               RNA.single.cell.type.specific.nTPM,
               regex(paste0("(?<=",cells,": )\\d+\\.\\d+"))),
             TPM = as.numeric(TPM)) %>% 
  dplyr::select(c(Ensembl, TPM)))
  })

GBM <- read.table(file.path(data, "cancer_glioma.tsv"),
                  header = T, sep = "\t") %>% 
  dplyr::mutate(
    GBM = str_extract(
      RNA.cancer.specific.FPKM,
      regex("(?<=glioma: )\\d+\\.\\d+")),
    GBM = as.numeric(GBM)) %>%
  dplyr::select(c(Ensembl, GBM))

cell_type_rna$GBMcells <- GBM

cell_type_rna <- merge.rec(cell_type_rna, by = "Ensembl", all = T,
                           suffixes=c("","")) %>%
  setNames(c("Ensembl", cell_types, "GBM_cells")) %>%
  dplyr::mutate(Ensembl = mapIds(org.Hs.eg.db, Ensembl,
                                 keytype = "ENSEMBL", column = "SYMBOL")) %>%
  dplyr::mutate_all(~replace(., is.na(.), 0)) %>% 
  dplyr::filter(Ensembl != 0) %>% 
  dplyr::group_by(Ensembl) %>%
  dplyr::summarise(across(everything(), sum)) %>% 
  tibble::column_to_rownames("Ensembl")

## Extract genes from the TxDb object
human.txdb <- makeTxDbFromGFF("A:/grch38/Homo_sapiens.GRCh38.gff", 
                              dataSource = "Ensembl", 
                              organism = "Homo sapiens")
# Extracting gene features from the TxDb object
human.genes <- exonsBy(human.txdb, by = "gene")
# Extracting gene lengths
exons <- reduce(human.genes)
geneLengths <- sum(width(exons))


TPM <- get_TPM(counts = assay(deseq$dds), effLen = geneLengths)
TPM <- TPM %>% as.data.frame() %>% 
  tibble::rownames_to_column("ensembl") %>% 
  dplyr::mutate(GeneSymbol = mapIds(org.Hs.eg.db, ensembl,
                                keytype = "ENSEMBL",column = "SYMBOL")) %>%
  dplyr::distinct(GeneSymbol, .keep_all = T) %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::select(!ensembl) %>% 
  tibble::column_to_rownames(., var="GeneSymbol")


decon <- deconvolute(m = as.matrix(TPM[,c(10,20)]),
                     sigMatrix = as.matrix(cell_type_rna))

# ------------------------------------------------------------------------------
# -                COMPARING PATIENT 593 AND PATIENT 673                       -
# ------------------------------------------------------------------------------
# 1.) Prepare the readcounts and metadata for the Deseq2 analysis
## split the readcounts and metadata for the two patients
readcounts.patients <- list("593" = readcounts[,grep("593", colnames(readcounts))],
                           "673" = readcounts[,grep("673", colnames(readcounts))])
coldata.patients <- coldata %>% 
  dplyr::mutate(condition = relevel(condition, ref = "physioxia.Tumour")) %>% 
  split.data.frame(.$patient, drop = T)

# 2.) Run the Deseq2 analysis
deseq.patients <- lapply(c("593","673"), function(x){
  make_deseq(matrix = readcounts.patients[[x]], 
             coldata = coldata.patients[[x]],
             design = "condition")
})
deseq.patients <- setNames(deseq.patients, c("593","673"))

# --------------------- Patient 593 --------------------------------------------
# Extract significant results for the 2D hypoxia vs 2D normoxia
results.593.2DN_vs_Tumor <- get_results(dds = deseq.patients$`593`$dds, # Deseq2 object
                                      sig_log2FC = 1.5, sig_pval = 0.05, # significance thresholds
                                      contrast=list("condition_normoxia.2D_vs_physioxia.Tumour"),
                                      name = "condition_normoxia.2D_vs_physioxia.Tumour")
# number of up- and downregulated genes: n = 476
table(results.593.2DN_vs_Tumor$sig_df$significance)
# Signif. down-regulated   Signif. up-regulated 
# 225                     251

# Extract significant results for the 3D organoids vs 2D normoxia
results.593.2DH_vs_Tumor <- get_results(dds = deseq.patients$`593`$dds, sig_log2FC = 1.5, sig_pval = 0.05,
                                     contrast=list("condition_hypoxia.2D_vs_physioxia.Tumour"),
                                     name = "condition_hypoxia.2D_vs_physioxia.Tumour")
# number of up- and downregulated genes: n = 435
table(results.593.2DH_vs_Tumor$sig_df$significance)
# Signif. down-regulated   Signif. up-regulated
# 193                      242

# Extract significant results for the Tumor vs 2D normoxia
results.593.3D_vs_Tumor <- get_results(dds = deseq.patients$`593`$dds, sig_log2FC = 1.5, sig_pval = 0.05,
                                        contrast=list("condition_physioxia.3D_vs_physioxia.Tumour"),
                                        name = "condition_physioxia.3D_vs_physioxia.Tumour")
# number of up- and downregulated genes: n = 496
table(results.593.3D_vs_Tumor$sig_df$significance)
# Signif. down-regulated   Signif. up-regulated
# 254                     242

## --------------------- Patient 673 -------------------------------------------
# Extract significant results for the 2D hypoxia vs 2D normoxia
results.673.2DN_vs_Tumor <- get_results(dds = deseq.patients$`673`$dds, 
                                      sig_log2FC = 1.5, sig_pval = 0.05,
                                      contrast=list("condition_normoxia.2D_vs_physioxia.Tumour"),
                                      name = "condition_normoxia.2D_vs_physioxia.Tumour")
# number of up- and downregulated genes: n = 702
table(results.673.2DN_vs_Tumor$sig_df$significance)
# Signif. down-regulated   Signif. up-regulated
# 353                     349

# Extract significant results for the 3D organoids vs 2D normoxia
results.673.2DH_vs_Tumor <- get_results(dds = deseq.patients$`673`$dds, 
                                        sig_log2FC = 1.5, sig_pval = 0.05,
                                     contrast=list("condition_hypoxia.2D_vs_physioxia.Tumour"),
                                     name = "condition_hypoxia.2D_vs_physioxia.Tumour")
# number of up- and downregulated genes: n = 543
table(results.673.2DH_vs_Tumor$sig_df$significance)
# Signif. down-regulated   Signif. up-regulated
# 248                      295

# Extract significant results for the Tumor vs 2D normoxia
results.673.3D_vs_Tumor <- get_results(dds = deseq.patients$`673`$dds, 
                                       sig_log2FC = 1.5, sig_pval = 0.05,
                                       contrast=list("condition_physioxia.3D_vs_physioxia.Tumour"),
                                       name = "condition_physioxia.3D_vs_physioxia.Tumour")
# number of up- and downregulated genes: n = 547
table(results.673.3D_vs_Tumor$sig_df$significance)
# Signif. down-regulated   Signif. up-regulated
# 312                      235

# Save results to excel files
sapply(c("df","sig_df"), function(x){
  # ----- Patient 593 -----
  # 2D normoxia vs Tumor
  write.csv(results.593.2DN_vs_Tumor[[x]],
            file.path(date, results_dir, paste0("593.2DN_vs_Tumor_",x,".csv")),
            row.names = T)
  # 2D hypoxia
  write.csv(results.593.2DH_vs_Tumor[[x]], 
            file.path(date, results_dir, paste0("593.2DH_vs_Tumor_",x,".csv")),
            row.names = T)
  # 3D organoids
  write.csv(results.593.3D_vs_Tumor[[x]], 
            file.path(date, results_dir, paste0("593.3D_vs_Tumor_",x,".csv")),
            row.names = T)
  # ----- Patient 673 -----
  # 2D normoxia vs Tumor
  write.csv(results.673.2DN_vs_Tumor[[x]],
            file.path(date, results_dir, paste0("673.2DN_vs_Tumor_",x,".csv")),
            row.names = T)
  # 2D hypoxia
  write.csv(results.673.2DH_vs_Tumor[[x]], 
            file.path(date, results_dir, paste0("673.2DH_vs_Tumor_",x,".csv")),
            row.names = T)
  # 3D organoids
  write.csv(results.673.3D_vs_Tumor[[x]], 
            file.path(date, results_dir, paste0("673.3D_vs_Tumor_",x,".csv")),
            row.names = T)
})

## ------------------------- Patient 593 ---------------------------------------
# KEGG pathways - ORA analysis
## 1. Get entrez IDs: gene set of interest and background
genes.593.2DN_vs_Tumor <- get_genelist(results.593.2DN_vs_Tumor$df, # 2D normoxia vs Tumor
                                     filter = results.593.2DN_vs_Tumor$sig_df$entrezID)
genes.593.2DH_vs_Tumor <- get_genelist(results.593.2DH_vs_Tumor$df, # 2D hypoxia vs Tumor
                                     filter = results.593.2DH_vs_Tumor$sig_df$entrezID)
genes.593.3D_vs_Tumor <- get_genelist(results.593.3D_vs_Tumor$df, # 3D organoids vs Tumor
                                     filter = results.593.3D_vs_Tumor$sig_df$entrezID)
## 2. Run the ORA analysis
KEGG.593.2DN_vs_Tumor <- list()
KEGG.593.2DN_vs_Tumor$df <- kegg_results(genes.593.2DN_vs_Tumor$interest,
                                         genes.593.2DN_vs_Tumor$background)

KEGG.593.2DH_vs_Tumor <- list()
KEGG.593.2DH_vs_Tumor$df <- kegg_results(genes.593.2DH_vs_Tumor$interest,
                                         genes.593.2DH_vs_Tumor$background)

KEGG.593.3D_vs_Tumor <- list()
KEGG.593.3D_vs_Tumor$df <- kegg_results(genes.593.3D_vs_Tumor$interest,
                                        genes.593.3D_vs_Tumor$background)
## 2.3 Save the results
openxlsx::write.xlsx(KEGG.593.2DN_vs_Tumor$df, 
                     file.path(date, results_dir, "593.2DN_vs_Tumor_KEGG.xlsx"))
openxlsx::write.xlsx(KEGG.593.2DH_vs_Tumor$df,
                     file.path(date, results_dir, "593.2DH_vs_Tumor_KEGG.xlsx"))
openxlsx::write.xlsx(KEGG.593.3D_vs_Tumor$df,
                     file.path(date, results_dir, "593.3D_vs_Tumor_KEGG.xlsx"))

## 1.2 GSEA GO analysis
GO.593.2DN_vs_Tumor <- list()
GO.593.2DN_vs_Tumor$df <- go_results(genes.593.2DN_vs_Tumor$interest,
                                   genes.593.2DN_vs_Tumor$background,
                                   type = 'ALL')
GO.593.2DH_vs_Tumor <- list()
GO.593.2DH_vs_Tumor$df <- go_results(genes.593.2DH_vs_Tumor$interest,
                                  genes.593.2DH_vs_Tumor$background,
                                  type = 'ALL')
GO.593.3D_vs_Tumor <- list()
GO.593.3D_vs_Tumor$df <- go_results(genes.593.3D_vs_Tumor$interest,
                                     genes.593.3D_vs_Tumor$background,
                                     type = 'ALL')
# Save the results
openxlsx::write.xlsx(GO.593.2DN_vs_Tumor$df,
                     file.path(date, results_dir, "593.2DN_vs_Tumor_GO.xlsx"))
openxlsx::write.xlsx(GO.593.2DH_vs_Tumor$df,
                     file.path(date, results_dir, "593.2DH_vs_Tumor_GO.xlsx"))
openxlsx::write.xlsx(GO.593.3D_vs_Tumor$df,
                     file.path(date, results_dir, "593.3D_vs_Tumor_GO.xlsx"))

## Combine the results
ORA.593.2DN_vs_Tumor <- bind_df(kegg = KEGG.593.2DN_vs_Tumor$df, go = GO.593.2DN_vs_Tumor) %>% 
  make_GObase(., results.593.2DN_vs_Tumor$sig_df)
ORA.593.2DH_vs_Tumor <- bind_df(kegg = KEGG.593.2DH_vs_Tumor$df, go = GO.593.2DH_vs_Tumor) %>% 
  make_GObase(., results.593.2DH_vs_Tumor$sig_df)
ORA.593.3D_vs_Tumor <- bind_df(kegg = KEGG.593.3D_vs_Tumor$df, go = GO.593.3D_vs_Tumor) %>% 
  make_GObase(., results.593.3D_vs_Tumor$sig_df)

## ------------------------- Patient 673 ---------------------------------------
# 2 KEGG pathways - ORA analysis
## 2.1 Get entrez IDs: gene set of interest and background
genes.673.2DN_vs_Tumor <- get_genelist(results.673.2DN_vs_Tumor$df, # 2D hypoxia vs 2D normoxia
                                     filter = results.673.2DN_vs_Tumor$sig_df$entrezID)
genes.673.2DH_vs_Tumor <- get_genelist(results.673.2DH_vs_Tumor$df, # 3D organoids vs 2D normoxia
                                    filter = results.673.2DH_vs_Tumor$sig_df$entrezID)
genes.673.3D_vs_Tumor <- get_genelist(results.673.3D_vs_Tumor$df, # Tumour vs 2D normoxia
                                       filter = results.673.3D_vs_Tumor$sig_df$entrezID)

## 2.2 Run the ORA analysis
KEGG.673.2DN_vs_Tumor <- list()
KEGG.673.2DN_vs_Tumor$df <- kegg_results(genes.673.2DN_vs_Tumor$interest,
                                       genes.673.2DN_vs_Tumor$background)

KEGG.673.2DH_vs_Tumor <- list()
KEGG.673.2DH_vs_Tumor$df <- kegg_results(genes.673.2DH_vs_Tumor$interest,
                                      genes.673.2DH_vs_Tumor$background)

KEGG.673.3D_vs_Tumor <- list()
KEGG.673.3D_vs_Tumor$df <- kegg_results(genes.673.3D_vs_Tumor$interest,
                                         genes.673.3D_vs_Tumor$background)

## 1.2 GSEA GO analysis
GO.673.2DN_vs_Tumor <- list()
GO.673.2DN_vs_Tumor$df <- go_results(genes.673.2DN_vs_Tumor$interest,
                                   genes.673.2DN_vs_Tumor$background,
                                   type = 'ALL')
GO.673.2DH_vs_Tumor <- list()
GO.673.2DH_vs_Tumor$df <- go_results(genes.673.2DH_vs_Tumor$interest,
                                  genes.673.2DH_vs_Tumor$background,
                                  type = 'ALL')
GO.673.3D_vs_Tumor <- list()
GO.673.3D_vs_Tumor$df <- go_results(genes.673.3D_vs_Tumor$interest,
                                     genes.673.3D_vs_Tumor$background,
                                     type = 'ALL')
# Save the results
openxlsx::write.xlsx(GO.673.2DN_vs_Tumor$df,
                     file.path(date, results_dir, "673.2DN_vs_Tumor_GO.xlsx"))
openxlsx::write.xlsx(GO.673.2DH_vs_Tumor$df,
                     file.path(date, results_dir, "673.2DH_vs_Tumor_GO.xlsx"))
openxlsx::write.xlsx(GO.673.3D_vs_Tumor$df,
                     file.path(date, results_dir, "673.3D_vs_Tumor_GO.xlsx"))



ORA.673.2DN_vs_Tumor <- bind_df(kegg = KEGG.673.2DN_vs_Tumor$df, go = GO.673.2DN_vs_Tumor) %>% 
  make_GObase(., results.673.2DN_vs_Tumor$sig_df)
ORA.673.2DH_vs_Tumor <- bind_df(kegg = KEGG.673.2DH_vs_Tumor$df, go = GO.673.2DH_vs_Tumor) %>% 
  make_GObase(., results.673.2DH_vs_Tumor$sig_df)
ORA.673.3D_vs_Tumor <- bind_df(kegg = KEGG.673.3D_vs_Tumor$df, go = GO.673.3D_vs_Tumor) %>% 
  make_GObase(., results.673.3D_vs_Tumor$sig_df)

####################################################
# 10.) Gene sets in interception                   #
####################################################
## ------------------------- Patient 593 ---------------------------------------
# Create Venn diagram of shared differentially expressed genes
simple_venn(results.593.2DN_vs_Tumor$sig_df, 
            results.593.2DH_vs_Tumor$sig_df, 
            results.593.3D_vs_Tumor$sig_df,
            c("2DN","2DH", "3D"),
            "venn_diagram_593")

venn.593 <- make_vennbase(results.593.2DN_vs_Tumor$sig_df, 
                          results.593.2DH_vs_Tumor$sig_df, 
                          results.593.3D_vs_Tumor$sig_df,
                          c("2DN","2DH", "3D"))

# get the membership of the genes in the regions
upset.593 <- make_upsetbase(venn.593$table, c("2DN","2DH", "3D"))
# get the centrum of the regions
arranged.593 <- arrange_venn(upset.593, c("2DN","2DH", "3D"),
                             extract_regions = T)

set.seed(123)
# randomize the position of the genes in the regions
xy.593 = rbind(
  calculateCircle(x = 0.00, y = 0.2886, r = 0.1, # 3-way intersection
                  noiseFun = function(x) (x + rnorm(1,0,0.1)), # add noise
                  steps = 207,randomDist = T, randomFun = rnorm), # 579 genes
  calculateEllipse(x = -0.65, y = -0.0866, a = 0.2, b = .2, angle = -240, 
                   noiseFun = function(x) (x + rnorm(1,0,0.1)),
                   steps = 53, randomDist = T, randomFun = rnorm),
  calculateEllipse(x = 0.65, y = -0.0866, a = 0.2, b = .2, angle = -120, 
                   noiseFun = function(x) (x + rnorm(1,0,0.1)),
                   steps = 16, randomDist = T, randomFun = rnorm),
  calculateEllipse(x = 0, y = 1.0392, a = 0.1, b = .2, angle = 0, 
                   noiseFun = function(x) (x + rnorm(1,0,0.1)),
                   steps = 160, randomDist = T, randomFun = rnorm),
  calculateCircle(x = 0.00, y = -1.2124, r = .3, 
                  noiseFun = function(x) (x + rnorm(1,0,0.2)),
                  steps = 220, randomDist = T, randomFun = rnorm),
  calculateCircle(x = 1.30, y = 1.0392, r = .3,  
                  noiseFun = function(x) (x + rnorm(1,0,0.2)),
                  steps = 52,randomDist = T, randomFun = rnorm),
  calculateCircle(x = -1.30, y = 1.0392, r = .3, 
                  noiseFun = function(x) (x + rnorm(1,0,0.2)),
                  steps = 56,randomDist = T, randomFun = rnorm)
)

upset.593 <- upset.593 %>%
  dplyr::mutate(region = dplyr::case_when(
    `2DN` & `2DH` & `3D` ~ "2DN-2DH-3D",
    `2DN` & !`2DH` & `3D` ~ "2DN-3D",
    !`2DN` & `2DH` & `3D` ~ "2DH-3D",
    `2DN` & `2DH` & !`3D` ~ "2DN-2DH",
    !`2DN` & !`2DH` & `3D` ~ "3D",
    !`2DN` & `2DH` & !`3D` ~ "2DH",
    `2DN` & !`2DH` & !`3D` ~ "2DN"
  ), # relabel the regions
  region = as.factor(region),
  x = xy.593[,1],
  y = xy.593[,2],
  ) %>% # add the x and y coordinates of the genes 
  # add their log2FoldChange values to the data frame
  dplyr::mutate(log2FoldChange = rowMeans(upset.593[,5:7], na.rm = T))

## ------------------------- Patient 673 ---------------------------------------
# Create Venn diagram of shared differentially expressed genes
simple_venn(results.673.2DN_vs_Tumor$sig_df, 
            results.673.2DH_vs_Tumor$sig_df, 
            results.673.3D_vs_Tumor$sig_df,
            c("2DN","2DH", "3D"),
            "venn_diagram_673")

venn.673 <- make_vennbase(results.673.2DN_vs_Tumor$sig_df, 
                          results.673.2DH_vs_Tumor$sig_df, 
                          results.673.3D_vs_Tumor$sig_df,
                          c("2DN","2DH", "3D"))

# get the membership of the genes in the regions
upset.673 <- make_upsetbase(venn.673$table, c("2DN", "2DH", "3D"))
# get the centrum of the regions
arranged.673 <- arrange_venn(upset.673, c("2DN", "2DH", "3D"),
                             extract_regions = T)

set.seed(123)
# randomize the position of the genes in the regions
xy.673 = rbind(
  calculateCircle(x = 0.00, y = 0.2886, r = 0.1, # 3-way intersection
                  noiseFun = function(x) (x + rnorm(1,0,0.1)), # add noise
                  steps = 265,randomDist = T, randomFun = rnorm), # 579 genes
  calculateEllipse(x = -0.65, y = -0.0866, a = 0.2, b = .2, angle = -240, 
                   noiseFun = function(x) (x + rnorm(1,0,0.1)),
                   steps = 61, randomDist = T, randomFun = rnorm),
  calculateEllipse(x = 0.65, y = -0.0866, a = 0.2, b = .2, angle = -120, 
                   noiseFun = function(x) (x + rnorm(1,0,0.1)),
                   steps = 8, randomDist = T, randomFun = rnorm),
  calculateEllipse(x = 0, y = 1.0392, a = 0.1, b = .2, angle = 0, 
                   noiseFun = function(x) (x + rnorm(1,0,0.1)),
                   steps = 252, randomDist = T, randomFun = rnorm),
  calculateCircle(x = 0.00, y = -1.2124, r = .3, 
                  noiseFun = function(x) (x + rnorm(1,0,0.2)),
                  steps = 213, randomDist = T, randomFun = rnorm),
  calculateCircle(x = 1.30, y = 1.0392, r = .3,  
                  noiseFun = function(x) (x + rnorm(1,0,0.2)),
                  steps = 18,randomDist = T, randomFun = rnorm),
  calculateCircle(x = -1.30, y = 1.0392, r = .3, 
                  noiseFun = function(x) (x + rnorm(1,0,0.2)),
                  steps = 124,randomDist = T, randomFun = rnorm)
)

upset.673 <- upset.673 %>%
  dplyr::mutate(region = dplyr::case_when(
    `2DN` & `2DH` & `3D` ~ "2DN-2DH-3D",
    `2DN` & !`2DH` & `3D` ~ "2DN-3D",
    !`2DN` & `2DH` & `3D` ~ "2DH-3D",
    `2DN` & `2DH` & !`3D` ~ "2DN-2DH",
    !`2DN` & !`2DH` & `3D` ~ "3D",
    !`2DN` & `2DH` & !`3D` ~ "2DH",
    `2DN` & !`2DH` & !`3D` ~ "2DN"
  ), # relabel the regions
  region = as.factor(region),
  x = xy.673[,1],
  y = xy.673[,2],
  ) %>% # add the x and y coordinates of the genes 
  # add their log2FoldChange values to the data frame
  dplyr::mutate(log2FoldChange = rowMeans(upset.673[,5:7], na.rm = T))

################################
# 12.) SURFME filter           #
################################

# 1.1 Prepare the data
surfme <- read.xlsx(file.path(data,"surfme_v2023.xlsx"))
surfme_categories <- c("GPI", "SinglePass", "MultiPass", "Cellular_Membrane",
                       "Extracellular_Domain", "GOCellSurf", "GOexternal","Surfy")
surfme_genes <- surfme %>% 
  dplyr::select(c(1,2,4, 6:13)) %>% 
  tidyr::separate_rows(Gene.names1, sep = " ") %>% 
  dplyr::rename("UniprotID" = Entry, 
                "Description" = Protein.names1,
                "geneID" = Gene.names1) %>% 
  dplyr::mutate(across(all_of(surfme_categories), ~ ifelse(is.na(.), F, T)))

surfme.593_all <- list(
  "2DN" = merge(results.593.2DN_vs_Tumor$df, surfme_genes, by = "geneID"),
  "2DH" = merge(results.593.2DH_vs_Tumor$df, surfme_genes, by = "geneID"),
  "3D" = merge(results.593.3D_vs_Tumor$df, surfme_genes, by = "geneID")
)

surfme.593_diff <- list(
  "2DN" = merge(results.593.2DN_vs_Tumor$sig_df, surfme_genes, by = "geneID"),
  "2DH" = merge(results.593.2DH_vs_Tumor$sig_df, surfme_genes, by = "geneID"),
  "3D" = merge(results.593.3D_vs_Tumor$sig_df, surfme_genes, by = "geneID")
)

write.xlsx(surfme.593_diff$`2DN`, file.path(date, results_dir,
                                                  "2DN_vs_Tumor_surfme.593.xlsx"))
write.xlsx(surfme.593_diff$`2DH`, file.path(date, results_dir,
                                                  "2DH_vs_Tumor_surfme.593.xlsx"))
write.xlsx(surfme.593_diff$`3D`, file.path(date, results_dir,
                                                 "3D_vs_Tumor_surfme.593.xlsx"))

surfme.593.upset_base <- list(
  data.frame("geneID" = surfme.593_diff$`2DN`$geneID,
             "2DN" = surfme.593_diff$`2DN`$log2FoldChange),
  data.frame("geneID" = surfme.593_diff$`2DH`$geneID,
             "2DH" = surfme.593_diff$`2DH`$log2FoldChange),
  data.frame("geneID" = surfme.593_diff$`3D`$geneID,
             "3D" = surfme.593_diff$`3D`$log2FoldChange))

surfme.593.upset_base <- merge.rec(surfme.593.upset_base, by="geneID", all=T, suffixes=c("", ""))
surfme.593.upset_base <- surfme.593.upset_base %>% 
  dplyr::distinct(., .keep_all = T) %>% 
  dplyr::mutate(across(all_of(c("X2DN","X2DH","X3D")), ~ ifelse(is.na(.), F, T))) %>% 
  dplyr::relocate(where(is.logical), .before = where(is.character)) %>% 
  setNames(c("2DN", "2DH", "3D", "geneID"))


surfme.593.upset_levels  <- as.factor(colnames(surfme.593.upset_base)[1:3])

###### Molecular functions
con <- GO.db::GO_dbconn()
go_term <- tbl(con, "go_term")
go_mf_parent <- tbl(con, "go_mf_parents")



MF_gene.593 <- AnnotationDbi::select(org.Hs.eg.db, surfme.593.upset_base$geneID, "GO","SYMBOL") %>%
  dplyr::filter(ONTOLOGY =="MF" & EVIDENCE != "ND") %>%
  dplyr::select(-c(EVIDENCE,ONTOLOGY)) %>%
  dplyr::rename("geneID" = SYMBOL) %>%
  distinct(geneID, .keep_all = T)

MF_gene.593 <- MF_gene.593 %>% 
  dplyr::mutate(
    AnnotationDbi::select(GO.db, GO, c("TERM"), "GOID")
  ) %>%
  dplyr::select(!GOID)


surfme.593.upset_MF <- merge(na.omit(surfme.593.upset_base), MF_gene.593, 
                        by = "geneID", all.x = T) %>%
  relocate(where(is.logical), .before = where(is.character)) %>% 
  dplyr::filter(!is.na(GO))

surfme.593.upset_MF <- na.omit(surfme.593.upset_MF) %>%
  dplyr::group_by(GO) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::ungroup()

other <- surfme.593.upset_MF %>% 
  dplyr::filter(count < 5) %>%
  nrow()
                
for (i in 1:10){
  print(paste0("Iteration: ", i))
  
  surfme.593.upset_MF <- surfme.593.upset_MF %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(
      GO = ifelse(count < 10, get_parent(GO, unique(surfme.593.upset_MF$GO)), GO)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(GO) %>%
    dplyr::mutate(count = n()) %>% 
    dplyr::ungroup()
  
  other <- c(other, surfme.593.upset_MF %>% 
               dplyr::filter(count < 10) %>%
               nrow())
}


MF_families.593 <- list(
  "protein binding" = c("protein binding","amide binding", "calcium ion binding", 
                           "carbohydrate derivative binding", "exogenous protein binding",
                           "ion binding", "nucleic acid binding"),
  "catalytic activity" = c("catalytic activity", "protein kinase activity",
                                           "endopeptidase inhibitor activity"),
  "transporter activity" = c("transporter activity", "chloride channel activity",
                                         "voltage-gated potassium channel activity",
                                         "sodium channel activity", "amino acid:sodium symporter activity"),
  "receptor activity" = c("transmembrane signaling receptor activity",
                                          "molecular transducer activity",
                                          "cargo receptor activity"),
  "structural activity" = c("structural molecule activity"),
  "other" = c("other")
)

surfme.593.upset_MF <- surfme.593.upset_MF %>%
  dplyr::mutate(
    AnnotationDbi::select(GO.db, GO, c("TERM"))
  ) %>% 
  dplyr::mutate(
    TERM = ifelse(count < 3, "other", TERM),
    TERM = factor(TERM, levels = unlist(MF_families)),
    TERM = relevel(TERM, "other")
  ) %>% 
  dplyr::select(!GOID)


MF_colors <- c(
  # Binding Activities
  "protein binding" = "#0000FF", "amide binding" = "#377eb8", 
  "calcium ion binding" = "#4ea3d1", "carbohydrate derivative binding" = "#6fb6e3",
  "exogenous protein binding" = "#91c9f5","ion binding"="#b3dcff", 
  "nucleic acid binding" = "lightblue",
  # Enzymatic and Catalytic Activities
  "catalytic activity" ="#008000", "protein kinase activity" = "green",
  "endopeptidase inhibitor activity" = "darkolivegreen2",
  #Transport and Channel Activities
  "transporter activity" = "darkred", "chloride channel activity" = "red", 
  "voltage-gated potassium channel activity" = "brown",
  "sodium channel activity" ="brown2", "amino acid:sodium symporter activity" = "#f24d50", 
  # Receptor and Signaling Activities
  "receptor activity" = "purple4",
  "transmembrane signaling receptor activity" = "magenta4",
  "molecular transducer activity" = "purple",
  "cargo receptor activity" = "violet", 
  # Structural Activities
  "structural molecule activity" = "orange","structural activity" = "orange",
  "other" = "#808080"
)

surfme.593.upset_MF <- surfme.593.upset_MF %>% 
  dplyr::mutate(super_family = sapply(TERM, map_to_super_family, super_families = MF_families),
                super_family = as.factor(super_family),
                super_family = relevel(super_family, "other"))

(upset_MF_plot.593 <- upset(data = surfme.593.upset_MF,
                       intersect = surfme.593.upset_levels, 
                       name='Conditions',
                       mode='distinct',
                       min_degree = 1,
                       max_degree = 3,
                       min_size = 5,
                       set_sizes = (
                         upset_set_size(position = "right") + 
                           geom_label(aes(label=..count..), stat='count', 
                                      position = position_stack(vjust = .5))
                       ),
                       base_annotations = list(
                         'Intersection size'=(
                           intersection_size(
                             counts = T,
                             bar_number_threshold = 1,
                             color = "black",
                             mapping=aes(
                               fill=TERM))
                           + scale_fill_manual(values = MF_colors)
                           + guides(fill = guide_legend(ncol = 1))
                           + theme(plot.background=element_rect(fill='#E5D3B3', color = "black",
                                                                size=0.5, linetype="solid"),
                                   legend.position = "right",
                                   legend.background = element_rect(fill="#E5D3B3", colour ="black",
                                                                    size=0.5, linetype="solid"),
                                   legend.spacing.y = unit(5,'mm'),
                                   legend.direction = "vertical",
                                   legend.key.size = unit(5,'mm'),
                                   legend.text = element_text(size=10))
                           + labs(y='# of genes', fill = "Molecular functions")
                         )
                       ),
                       annotations =list(
                         'Percentage of Group'=list(
                           aes=aes(x=intersection, group = super_family,
                                   fill=super_family),
                           geom=list(
                             geom_bar(stat='count', position='fill', 
                                      show.legend = F),
                             geom_text(
                               aes(
                                 label=!!aes_percentage(relative_to='intersection'),
                                 group = super_family
                               ),
                               stat='count',
                               position=position_fill(vjust = .5)
                             ),
                             scale_y_continuous(labels=scales::percent_format()),
                             scale_fill_manual(values = MF_colors)
                           )
                         )
                       ),
                       themes=upset_modify_themes(
                         list(
                           'intersections_matrix'=theme(text=element_text(size=16)),
                           'overall_sizes'=theme(text=element_text(size=16))
                         )),
                       stripes=c('#E5D3B3', 'white')))

ggsave(file.path(date, plots_dir, "png","surfme.593_MF_upset.png"), "png", 
       plot = upset_MF_plot, width = 16, height = 16, dpi = 300)
ggsave(file.path(date, plots_dir, "svg","surfme.593_MF_upset.svg"), "svg", 
       plot = upset_MF_plot, width = 20, height = 16, units = 'in')


# ------------------------ Patient 673 ---------------------------------------
surfme.673_all <- list(
  "2DN" = merge(results.673.2DN_vs_Tumor$df, surfme_genes, by = "geneID"),
  "2DH" = merge(results.673.2DH_vs_Tumor$df, surfme_genes, by = "geneID"),
  "3D" = merge(results.673.3D_vs_Tumor$df, surfme_genes, by = "geneID")
)

surfme.673_diff <- list(
  "2DN" = merge(results.673.2DN_vs_Tumor$sig_df, surfme_genes, by = "geneID"),
  "2DH" = merge(results.673.2DH_vs_Tumor$sig_df, surfme_genes, by = "geneID"),
  "3D" = merge(results.673.3D_vs_Tumor$sig_df, surfme_genes, by = "geneID")
)

write.xlsx(surfme.673_diff$`2DN`, file.path(date, results_dir,
                                            "2DN_vs_Tumor_surfme.673.xlsx"))
write.xlsx(surfme.673_diff$`2DH`, file.path(date, results_dir,
                                            "2DH_vs_Tumor_surfme.673.xlsx"))
write.xlsx(surfme.673_diff$`3D`, file.path(date, results_dir,
                                           "3D_vs_Tumor_surfme.673.xlsx"))

surfme.673.upset_base <- list(
  data.frame("geneID" = surfme.673_diff$`2DN`$geneID,
             "2DN" = surfme.673_diff$`2DN`$log2FoldChange),
  data.frame("geneID" = surfme.673_diff$`2DH`$geneID,
             "2DH" = surfme.673_diff$`2DH`$log2FoldChange),
  data.frame("geneID" = surfme.673_diff$`3D`$geneID,
             "3D" = surfme.673_diff$`3D`$log2FoldChange))

surfme.673.upset_base <- merge.rec(surfme.673.upset_base, by="geneID", all=T, suffixes=c("", ""))
surfme.673.upset_base <- surfme.673.upset_base %>% 
  dplyr::distinct(., .keep_all = T) %>% 
  dplyr::mutate(across(all_of(c("X2DN","X2DH","X3D")), ~ ifelse(is.na(.), F, T))) %>% 
  dplyr::relocate(where(is.logical), .before = where(is.character)) %>% 
  setNames(c("2DN", "2DH", "3D", "geneID"))


surfme.673.upset_levels  <- as.factor(colnames(surfme.673.upset_base)[1:3])

MF_gene.673 <- AnnotationDbi::select(org.Hs.eg.db, surfme.673.upset_base$geneID, "GO","SYMBOL") %>%
  dplyr::filter(ONTOLOGY =="MF" & EVIDENCE != "ND") %>%
  dplyr::select(-c(EVIDENCE,ONTOLOGY)) %>%
  dplyr::rename("geneID" = SYMBOL) %>%
  distinct(geneID, .keep_all = T)

MF_gene.673 <- MF_gene.673 %>% 
  dplyr::mutate(
    AnnotationDbi::select(GO.db, GO, c("TERM"), "GOID")
  ) %>%
  dplyr::select(!GOID)


surfme.673.upset_MF <- merge(na.omit(surfme.673.upset_base), MF_gene.673, 
                             by = "geneID", all.x = T) %>%
  relocate(where(is.logical), .before = where(is.character)) %>% 
  dplyr::filter(!is.na(GO))

surfme.673.upset_MF <- na.omit(surfme.673.upset_MF) %>%
  dplyr::group_by(GO) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::ungroup()

other <- surfme.673.upset_MF %>% 
  dplyr::filter(count < 5) %>%
  nrow()

for (i in 1:10){
  print(paste0("Iteration: ", i))
  
  surfme.673.upset_MF <- surfme.673.upset_MF %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(
      GO = ifelse(count < 5, get_parent(GO, unique(surfme.673.upset_MF$GO)), GO)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(GO) %>%
    dplyr::mutate(count = n()) %>% 
    dplyr::ungroup()
  
  other <- c(other, surfme.673.upset_MF %>% 
               dplyr::filter(count < 5) %>%
               nrow())
}


MF_families.673 <- list(
  "protein binding" = c("protein binding","amide binding", "calcium ion binding", 
                        "exogenous protein binding"),
  "catalytic activity" = c("catalytic activity", "protein kinase activity",
                           "endopeptidase inhibitor activity", "molecular function regulator activity"),
  "transporter activity" = c("transporter activity", "chloride channel activity",
                             "voltage-gated potassium channel activity"),
  "receptor activity" = c("transmembrane receptor protein tyrosine kinase activity",
                          "molecular transducer activity", "G protein-coupled receptor activity"),
  "structural activity" = c("structural molecule activity"),
  "other" = c("other")
)

surfme.673.upset_MF <- surfme.673.upset_MF %>%
  dplyr::mutate(
    AnnotationDbi::select(GO.db, GO, c("TERM"))
  ) %>% 
  dplyr::group_by(GO) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(
    TERM = ifelse(count < 3, "other", TERM),
    TERM = factor(TERM, levels = unlist(MF_families.673)),
    TERM = relevel(TERM, "other")
  ) %>% 
  dplyr::select(!GOID)


MF_colors <- c(
  # Binding Activities
  "protein binding" = "#0000FF", "amide binding" = "#377eb8", 
  "calcium ion binding" = "#4ea3d1", "exogenous protein binding" = "#91c9f5",
  # Enzymatic and Catalytic Activities
  "catalytic activity" ="#008000", "protein kinase activity" = "green",
  "endopeptidase inhibitor activity" = "darkolivegreen2", 
  "molecular function regulator activity" = "lightgreen",
  #Transport and Channel Activities
  "transporter activity" = "darkred", "chloride channel activity" = "red", 
  "voltage-gated potassium channel activity" = "brown2", 
  # Receptor and Signaling Activities
  "receptor activity" = "purple4",
  "transmembrane receptor protein tyrosine kinase activity" = "magenta4",
  "molecular transducer activity" = "purple",
  "G protein-coupled receptor activity" = "violet", 
  # Structural Activities
  "structural molecule activity" = "orange","structural activity" = "orange",
  "other" = "#808080"
)

surfme.673.upset_MF <- surfme.673.upset_MF %>% 
  dplyr::mutate(super_family = sapply(TERM, map_to_super_family, super_families = MF_families.673),
                super_family = as.factor(super_family),
                super_family = relevel(super_family, "other"))

(upset_MF_plot.673 <- upset(data = surfme.673.upset_MF,
                            intersect = surfme.673.upset_levels, 
                            name='Conditions',
                            mode='distinct',
                            min_degree = 1,
                            max_degree = 3,
                            min_size = 1,
                            set_sizes = (
                              upset_set_size(position = "right") + 
                                geom_label(aes(label=..count..), stat='count', 
                                           position = position_stack(vjust = .5))
                            ),
                            base_annotations = list(
                              'Intersection size'=(
                                intersection_size(
                                  counts = T,
                                  bar_number_threshold = 1,
                                  color = "black",
                                  mapping=aes(
                                    fill=TERM))
                                + scale_fill_manual(values = MF_colors)
                                + guides(fill = guide_legend(ncol = 1))
                                + theme(plot.background=element_rect(fill='#E5D3B3', color = "black",
                                                                     size=0.5, linetype="solid"),
                                        legend.position = "right",
                                        legend.background = element_rect(fill="#E5D3B3", colour ="black",
                                                                         size=0.5, linetype="solid"),
                                        legend.spacing.y = unit(5,'mm'),
                                        legend.direction = "vertical",
                                        legend.key.size = unit(5,'mm'),
                                        legend.text = element_text(size=10))
                                + labs(y='# of genes', fill = "Molecular functions")
                              )
                            ),
                            annotations =list(
                              'Percentage of Group'=list(
                                aes=aes(x=intersection, group = super_family,
                                        fill=super_family),
                                geom=list(
                                  geom_bar(stat='count', position='fill', 
                                           show.legend = F),
                                  geom_text(
                                    aes(
                                      label=!!aes_percentage(relative_to='intersection'),
                                      group = super_family
                                    ),
                                    stat='count',
                                    position=position_fill(vjust = .5)
                                  ),
                                  scale_y_continuous(labels=scales::percent_format()),
                                  scale_fill_manual(values = MF_colors)
                                )
                              )
                            ),
                            themes=upset_modify_themes(
                              list(
                                'intersections_matrix'=theme(text=element_text(size=16)),
                                'overall_sizes'=theme(text=element_text(size=16))
                              )),
                            stripes=c('#E5D3B3', 'white')))

ggsave(file.path(date, plots_dir, "png","surfme.673_MF_upset.png"), "png", 
       plot = upset_MF_plot, width = 16, height = 16, dpi = 300)
ggsave(file.path(date, plots_dir, "svg","surfme.673_MF_upset.svg"), "svg", 
       plot = upset_MF_plot, width = 20, height = 16, units = 'in')


surfme_unchanged <- setdiff(c(surfme_all$`2DN`$geneID, surfme_all$`2DH`$geneID, surfme_all$`3D`$geneID),
                            c(surfme_diff$`2DN`$geneID, surfme_diff$`2DH`$geneID, surfme_diff$`3D`$geneID))
full_expr <- list(
  "2DN" = results.2DN_vs_Tumor$df %>% select("geneID", "baseMean", "log2FoldChange", "significance"),
  "2DH" = results.2DH_vs_Tumor$df %>% select("geneID", "baseMean", "log2FoldChange", "significance"),
  "3D" = results.3D_vs_Tumor$df %>% select("geneID", "baseMean", "log2FoldChange", "significance")
) %>% merge.rec(by = "geneID", all = T, suffixes = c("", ""))

surfme_unchanged_df <- subset.data.frame(full_expr, geneID %in% surfme_unchanged) %>%
  merge(., surfme_genes, by = "geneID", all.x = T)

surfme_unchanged_df <- surfme_unchanged_df %>% 
  dplyr::filter(if_all(starts_with("significance"), ~ . == "NS")) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(meanExpr = mean(c_across(starts_with("baseMean")), na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(geneID, Description, UniprotID, meanExpr, 
                log2FoldChange, log2FoldChange.1, log2FoldChange.2,
                all_of(surfme_categories))

write.xlsx(surfme_unchanged_df, file.path(date, results_dir, "unchanged_surfme.xlsx"))

