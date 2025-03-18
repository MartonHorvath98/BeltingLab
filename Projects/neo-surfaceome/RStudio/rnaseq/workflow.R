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
  mrna.path <- file.path("../../results/04_counts/")
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
    # clean up the column names
    setNames(c("Geneid", gsub(".counts.txt","",file.names))) %>% 
    # Clean up gene names
    dplyr::mutate(Geneid = gsub("\\..*","",Geneid)) %>%
    # set the gene ID as row names
    dplyr::distinct(Geneid, .keep_all = T) %>%
    tibble::column_to_rownames("Geneid") 
  
  write.table(mrna.counts, paste(data,"readcounts.csv", sep = "/"), 
              sep =",", na = "NA", dec = ".", row.names = T, col.names = T)
  
  # Load count read counts matrix
  readcounts <- read.csv(file = file.path(data,"readcounts.csv"),
                         sep = ",", header = T, na.strings = NA, row.names = 1)
  readcounts <- as.matrix(readcounts)
  
}
# Calculate TPM values for each sample
if (exists("TPM") == F){
  # Create a TxDb object from the GTF reference annotation
  human.txdb <- txdbmaker::makeTxDbFromGFF("../../refs/ref_annot.gtf", format = "gtf",                              format = "gtf",
                                dataSource = "Ensembl", 
                                organism = "Homo sapiens")
  # Extract gene features from the TxDb object
  human.genes <- exonsBy(human.txdb, by = "gene")
  # Extract the gene lengths
  exons <- reduce(human.genes)
  geneLengths <- sum(width(exons))
  # Calculate TPM values
  TPM <- get_TPM(counts = readcounts, effLen = geneLengths)
  TPM <- TPM[which(rowSums(cpm(readcounts, normalized = T) >1 ) >= 3),]
  ### Quantile normalize TPM value
  TPM <- normalize.quantiles(TPM, keep.names = T)
  row.names(TPM) <- row.names(readcounts)
  
  write.table(TPM, paste(data,"TPM.csv", sep = "/"), 
              sep =",", na = "NA", dec = ".", row.names = T, col.names = T)
  # Load count read counts matrix
  TPM <- read.csv(file = file.path(data,"TPM.csv"),
                         sep = ",", header = T, na.strings = NA, row.names = 1)
  TPM <- as.matrix(TPM)
}
#clean up the workspace
rm(list = c("mrna.reads","mrna.counts","file.names", "mrna.path", "mrna.files"))

# ---------------------------------------------------------------------------- #
# -            RUN DIFFERENTIAL EXPRESSION ANALYSIS                          - #
# ---------------------------------------------------------------------------- #

if (exists("deseq") == F) {
  # 1.) Prepare metadata table based on the sample names
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
  
  # 3.) split the readcounts and metadata for the two patients
  readcounts.patients <- list("593" = readcounts[,grep("593", colnames(readcounts))],
                              "673" = readcounts[,grep("673", colnames(readcounts))])
  coldata.patients <- coldata %>% 
    dplyr::mutate(condition = relevel(condition, ref = "normoxia.2D")) %>% 
    split.data.frame(.$patient, drop = T)
  
  # 2.) Run the Deseq2 analysis
  deseq.patients <- lapply(c("593","673"), function(x){
    make_deseq(matrix = readcounts.patients[[x]], 
               coldata = coldata.patients[[x]],
               design = "condition")
  })
  deseq.patients <- setNames(deseq.patients, c("593","673"))
}
# Clean up the workspace
rm(list = c("samples","human.txdb","human.genes","exons","geneLengths"))
# ------------------------------------------------------------------------------
# -                COMPARING PATIENT 593 AND PATIENT 673                       -
# ------------------------------------------------------------------------------

if (exists("results.593") == F | exists("results.673")) {
  # --------------------- Patient 593 --------------------------------------------
  # Extract significant results for the 2D hypoxia vs 2D normoxia
  results.593 <- list()
  results.593[["2DN_vs_2DH"]] <- get_results(dds = deseq.patients$`593`$dds, # Deseq2 object
                                          sig_log2FC = 1.5, sig_pval = 0.05, # significance thresholds
                                          contrast=list("condition_hypoxia.2D_vs_normoxia.2D"),
                                          name = "condition_hypoxia.2D_vs_normoxia.2D")
  # number of up- and downregulated genes: n = 113
  table(results.593$`2DN_vs_2DH`$sig_df$significance)
  # Signif. down-regulated   Signif. up-regulated 
  # 6                        107
  
  # Extract significant results for the 3D organoids vs 2D normoxia
  results.593[["2DN_vs_3D"]] <- get_results(dds = deseq.patients$`593`$dds, sig_log2FC = 1.5, sig_pval = 0.05,
                                          contrast=list("condition_physioxia.3D_vs_normoxia.2D" ),
                                          name = "condition_physioxia.3D_vs_normoxia.2D" )
  # number of up- and downregulated genes: n = 258
  table(results.593$`2DN_vs_3D`$sig_df$significance)
  # Signif. down-regulated   Signif. up-regulated
  # 166                      92
  
  # Extract significant results for the Tumor vs 2D normoxia
  results.593[["2DN_vs_Tumor"]] <- get_results(dds = deseq.patients$`593`$dds, sig_log2FC = 1.5, sig_pval = 0.05,
                                         contrast=list("condition_physioxia.Tumour_vs_normoxia.2D"),
                                         name = "condition_physioxia.Tumour_vs_normoxia.2D")
  # number of up- and downregulated genes: n = 466
  table(results.593$`2DN_vs_Tumor`$sig_df$significance)
  # Signif. down-regulated   Signif. up-regulated
  # 245                      221  
  
  ## --------------------- Patient 673 -------------------------------------------
  results.673 <- list()
  # Extract significant results for the 2D hypoxia vs 2D normoxia
  results.673[["2DN_vs_2DH"]] <- get_results(dds = deseq.patients$`673`$dds, 
                                          sig_log2FC = 1.5, sig_pval = 0.05,
                                          contrast=list("condition_hypoxia.2D_vs_normoxia.2D"),
                                          name = "condition_hypoxia.2D_vs_normoxia.2D")
  # number of up- and downregulated genes: n = 86
  table(results.673$`2DN_vs_2DH`$sig_df$significance)
  # Signif. down-regulated   Signif. up-regulated
  # 1                        85
  
  # Extract significant results for the 3D organoids vs 2D normoxia
  results.673[["2DN_vs_3D"]] <- get_results(dds = deseq.patients$`673`$dds, 
                                          sig_log2FC = 1.5, sig_pval = 0.05,
                                          contrast=list("condition_physioxia.3D_vs_normoxia.2D"),
                                          name = "condition_physioxia.3D_vs_normoxia.2D")
  # number of up- and downregulated genes: n = 792
  table(results.673$`2DN_vs_3D`$sig_df$significance)
  # Signif. down-regulated   Signif. up-regulated
  # 473                      319
  
  # Extract significant results for the Tumor vs 2D normoxia
  results.673[["2DN_vs_Tumor"]] <- get_results(dds = deseq.patients$`673`$dds, 
                                         sig_log2FC = 1.5, sig_pval = 0.05,
                                         contrast=list("condition_physioxia.Tumour_vs_normoxia.2D"),
                                         name = "condition_physioxia.Tumour_vs_normoxia.2D")
  # number of up- and downregulated genes: n = 768
  table(results.673$`2DN_vs_Tumor`$sig_df$significance)
  # Signif. down-regulated   Signif. up-regulated
  # 393                      375
  
  # Save results to excel files
  sapply(c("df","sig_df"), function(x){
    # ----- Patient 593 -----
    # 2D hypoxia vs nromoxia
    write.xlsx(results.593$`2DN_vs_2DH`[[x]], 
               file = file.path(date, results_dir, "Patien593_2DN_vs_2DH.xlsx"), 
               sheetName = x, row.names = F)
    # 3D organoids vs 2D normoxia
    write.xlsx(results.593$`2DN_vs_3D`[[x]], 
               file = file.path(date, results_dir, "Patient593_2DN_vs_3D.xlsx"), 
               sheetName = x, row.names = F)
    # Tumor vs 2D normoxia
    write.xlsx(results.593$`2DN_vs_Tumor`[[x]], 
               file = file.path(date, results_dir, "Patient593_2DN_vs_Tumor.xlsx"), 
               sheetName = x, row.names = F)
    # ----- Patient 673 -----
    # 2D hypoxia vs nromoxia
    write.xlsx(results.673$`2DN_vs_2DH`[[x]], 
               file = file.path(date, results_dir, "Patient673_2DN_vs_2DH.xlsx"), 
               sheetName = x, row.names = F)
    # 3D organoids vs 2D normoxia
    write.xlsx(results.673$`2DN_vs_3D`[[x]], 
               file = file.path(date, results_dir, "Patient673_2DN_vs_3D.xlsx"), 
               sheetName = x, row.names = F)
    # Tumor vs 2D normoxia
    write.xlsx(results.673$`2DN_vs_Tumor`[[x]], 
               file = file.path(date, results_dir, "Patient673_2DN_vs_Tumor.xlsx"), 
               sheetName = x, row.names = F)
  })
}
# ---------------------------------------------------------------------------- #
# -       Functional analyses: ORA and GSEA with pathways and terms          - #
# ---------------------------------------------------------------------------- #
if (exists("GO.673") == F) {
  comparisons <- c("2DN_vs_2DH","2DN_vs_3D","2DN_vs_Tumor")
  # 1.) Extract entrez IDs for gene set of interest and background
  ### ------- Patient 593 -------
  gene_lists.593 <- list()
  gene_lists.593 <- lapply(
    list(results.593$`2DN_vs_2DH`$df,results.593$`2DN_vs_3D`$df,results.593$`2DN_vs_Tumor`$df),
    function(x){
      get_genelist(.df = x, 
                   .filter = x[["significance"]] %in% c("Signif. up-regulated", 
                                                      "Signif. down-regulated"))
  })
  names(gene_lists.593) <- comparisons
  
  ### ------- Patient 673 --------
  gene_lists.673 <- list()
  gene_lists.673 <- lapply(
    list(results.673$`2DN_vs_2DH`$df,results.673$`2DN_vs_3D`$df,results.673$`2DN_vs_Tumor`$df),
    function(x){
      get_genelist(.df = x, 
                   .filter = x[["significance"]] %in% c("Signif. up-regulated", 
                                                      "Signif. down-regulated"))
  })
  names(gene_lists.673) <- comparisons
  
  # 2.) Run the overrepresentation analysis (ORA)
  ### ------- Patient 593 -------
  ORA.593 <- list()
  ORA.593 <- lapply(gene_lists.593,
                     function (x){
                       run_ora(.interest = x[["interest"]],
                               .background = x[["background"]],
                               .pathways = pathways)})
  # Extract results
  for (i in names(ORA.593)){
    ORA.593[[i]] <- c(ORA.593[[i]],
                      extract_ora_results(.ora = ORA.593[[i]]$ora,
                                          .db = pathways,
                                          .expr = results.593[[i]]$df))
  }
  # Save the results
  sapply(names(ORA.593), function(x){
    openxlsx::write.xlsx(ORA.593[[x]][c(2:3)], 
                         file.path(date, results_dir, paste0("Patient593_", x, "_pathways_ORA.xlsx")))
  })
  ### ------- Patient 673 --------
  ORA.673 <- list()
  ORA.673 <- lapply(gene_lists.673,
                     function (x){
                       run_ora(.interest = x[["interest"]],
                               .background = x[["background"]],
                               .pathways = pathways)})
  # Extract results
  for (i in names(ORA.673)){
    ORA.673[[i]] <- c(ORA.673[[i]],
                      extract_ora_results(.ora = ORA.673[[i]]$ora,
                                          .db = pathways,
                                          .expr = results.673[[i]]$df))
  }
  
  # Save the results
  sapply(names(ORA.673), function(x){
    openxlsx::write.xlsx(ORA.673[[x]][c(2:3)], 
                         file.path(date, results_dir, paste0("Patient673_", x, "_pathways_ORA.xlsx")))
  })
  
  
  # 3.) Combine and cluster ORA results 
  ORA.593.compare <- list(ORA.593$`2DN_vs_2DH`$sig_df[,c("Database","ID","Description", "zscore", "p.adjust","Count", "geneID")],
                          ORA.593$`2DN_vs_3D`$sig_df[,c("Database","ID","Description", "zscore", "p.adjust","Count", "geneID")],
                          ORA.593$`2DN_vs_Tumor`$sig_df[,c("Database","ID","Description", "zscore", "p.adjust","Count", "geneID")]) %>%
    merge.rec(., by = c("Database","ID","Description"),
              suffixes = c("",""), all = T) %>% 
    setNames(., c("Database","ID","Description",
                  paste(rep(c("zscore","padj","Count","core_enrichment"),3),
                        rep(c("2D","3D","Tumor"),times=c(4,4,4)),
                        sep = ".")))
  ORA.673.compare <- list(ORA.673$`2DN_vs_2DH`$sig_df[,c("Database","ID","Description", "zscore", "p.adjust","Count", "geneID")],
                          ORA.673$`2DN_vs_3D`$sig_df[,c("Database","ID","Description", "zscore", "p.adjust","Count", "geneID")],
                          ORA.673$`2DN_vs_Tumor`$sig_df[,c("Database","ID","Description", "zscore", "p.adjust","Count", "geneID")]) %>%
    merge.rec(., by = c("Database","ID","Description"),
              suffixes = c("",""), all = T) %>% 
    setNames(., c("Database","ID","Description",
                  paste(rep(c("zscore","padj","Count","core_enrichment"),3),
                        rep(c("2D","3D","Tumor"),times=c(4,4,4)),
                        sep = ".")))
  
  ORA.593.compare <- cluster_enrichment(ORA.593.compare, pathways, "core_enrichment")
  ORA.673.compare <- cluster_enrichment(ORA.673.compare, pathways, "core_enrichment")
  
  # 4.) Run the gene set enrichment analysis (GSEA)
  ### ------- Patient 593 -------
  GO.593 <- list()
  GO.593 <- lapply(gene_lists.593, function(x){
    run_gsea(.geneset = x[["background"]], .terms = terms)})
  for (i in names(GO.593)){
    GO.593[[i]] <- c(GO.593[[i]],
                      extract_gsea_results(.gsea = GO.593[[i]]$gsea,
                                          .db = terms,
                                          .expr = results.593[[i]]$df))
  }
  # Save the results
  sapply(names(GO.593), function(x){
    openxlsx::write.xlsx(GO.593[[x]][c(2:3)], 
                         file.path(date, results_dir, paste0("Patient593_", x, "_go&hallmark_GSEA.xlsx")))
  })
  ### ------- Patient 673 --------
  GO.673 <- list()
  GO.673 <- lapply(gene_lists.673, function(x){
    run_gsea(.geneset = x[["background"]], .terms = terms)})
  for (i in names(GO.673)){
    GO.673[[i]] <- c(GO.673[[i]],
                     extract_gsea_results(.gsea = GO.673[[i]]$gsea,
                                          .db = terms,
                                          .expr = results.673[[i]]$df))
  }
  # Save the results
  sapply(names(GO.673), function(x){
    openxlsx::write.xlsx(GO.673[[x]][c(2:3)], 
                         file.path(date, results_dir, paste0("Patient673_", x, "_go&hallmark_GSEA.xlsx")))
  })
  
  # 5.) Combine and cluster GSEA results 
  GO.593.compare <- list(GO.593$`2DN_vs_2DH`$sig_df[,c("Database","ID", "NES", "p.adjust","zscore", "core_enrichment")] %>% 
                           tibble::rownames_to_column("Name"),
                         GO.593$`2DN_vs_3D`$sig_df[,c("Database","ID", "NES", "p.adjust","zscore", "core_enrichment")] %>% 
                           tibble::rownames_to_column("Name"),
                         GO.593$`2DN_vs_Tumor`$sig_df[,c("Database","ID", "NES", "p.adjust","zscore", "core_enrichment")] %>% 
                           tibble::rownames_to_column("Name")) %>%
    merge.rec(., by = c("Database","ID","Name"),
              suffixes = c("",""), all = T) %>% 
    setNames(., c("Database","ID","Description",
                  paste(rep(c("NES","padj","zscore","core_enrichment"),3),
                        rep(c("2D","3D","Tumor"),times=c(4,4,4)),
                        sep = ".")))
  GO.673.compare <- list(GO.673$`2DN_vs_2DH`$sig_df[,c("Database","ID", "NES", "p.adjust","zscore", "core_enrichment")] %>% 
                            tibble::rownames_to_column("Name"),
                          GO.673$`2DN_vs_3D`$sig_df[,c("Database","ID", "NES", "p.adjust","zscore", "core_enrichment")] %>% 
                            tibble::rownames_to_column("Name"),
                          GO.673$`2DN_vs_Tumor`$sig_df[,c("Database","ID", "NES", "p.adjust","zscore", "core_enrichment")] %>% 
                            tibble::rownames_to_column("Name")) %>%
    merge.rec(., by = c("Database","ID","Name"),
              suffixes = c("",""), all = T) %>% 
    setNames(., c("Database","ID","Description",
                  paste(rep(c("NES","padj","zscore","core_enrichment"),3),
                        rep(c("2D","3D","Tumor"),times=c(4,4,4)),
                        sep = ".")))
  
  
  GO.593.compare <- cluster_enrichment(GO.593.compare, terms, "core_enrichment")
  GO.673.compare <- cluster_enrichment(GO.673.compare, terms, "core_enrichment")
  
}


# ---------------------------------------------------------------------------- #
# -             Venn diagrams from gene sets in interception                 - #
# ---------------------------------------------------------------------------- #
if (exists("venn.593") == F) {
  # 1.) Create Venn diagram of shared differential expressed genes
  
  ### ------------------------- Patient 593 ------------------------------------
  venn.593 <- make_vennbase(results.593$`2DN_vs_2DH`$sig_df, 
                            results.593$`2DN_vs_3D`$sig_df, 
                            results.593$`2DN_vs_Tumor`$sig_df,
                            c("2DH","3D", "Tumor"))
  # get the membership of the genes in the regions
  upset.593 <- make_upsetbase(venn.593$table, c("2DH","3D", "Tumor"))
  # get the centrum of the regions
  arranged.593 <- arrange_venn(upset.593, c("2DH","3D", "Tumor"),
                               extract_regions = T)
  # randomize the position of the genes in the regions
  set.seed(42)
  xy.593 = rbind(
    calculateCircle(x = 0.00, y = 0.2886, r = 0.1, # 3-way intersection
                    noiseFun = function(x) (x + rnorm(1,0,0.1)), # add noise
                    steps = 2,randomDist = T, randomFun = rnorm), # 2 genes
    calculateEllipse(x = -0.65, y = -0.0866, a = 0.2, b = .2, angle = -240, 
                     noiseFun = function(x) (x + rnorm(1,0,0.1)),
                     steps = 22, randomDist = T, randomFun = rnorm),
    calculateEllipse(x = 0.65, y = -0.0866, a = 0.2, b = .2, angle = -120, 
                     noiseFun = function(x) (x + rnorm(1,0,0.1)),
                     steps = 76, randomDist = T, randomFun = rnorm),
    calculateEllipse(x = 0, y = 1.0392, a = 0.1, b = .2, angle = 0, 
                     noiseFun = function(x) (x + rnorm(1,0,0.1)),
                     steps = 2, randomDist = T, randomFun = rnorm),
    calculateCircle(x = 0.00, y = -1.2124, r = .3, 
                    noiseFun = function(x) (x + rnorm(1,0,0.2)),
                    steps = 366, randomDist = T, randomFun = rnorm),
    calculateCircle(x = 1.30, y = 1.0392, r = .3,  
                    noiseFun = function(x) (x + rnorm(1,0,0.2)),
                    steps = 178,randomDist = T, randomFun = rnorm),
    calculateCircle(x = -1.30, y = 1.0392, r = .3, 
                    noiseFun = function(x) (x + rnorm(1,0,0.2)),
                    steps = 87,randomDist = T, randomFun = rnorm)
  )
  # create data table for plotting
  upset.593 <- upset.593 %>%
    dplyr::mutate(region = dplyr::case_when(
      `2DH` & `3D` & `Tumor` ~ "2DH-3D-Tumor",
      `2DH` & !`3D` & `Tumor` ~ "2DH-Tumor",
      !`2DH` & `3D` & `Tumor` ~ "3D-Tumor",
      `2DH` & `3D` & !`Tumor` ~ "2DH-3D",
      !`2DH` & !`3D` & `Tumor` ~ "Tumor",
      !`2DH` & `3D` & !`Tumor` ~ "3D",
      `2DH` & !`3D` & !`Tumor` ~ "2DH"
    ), # relabel the regions
    region = as.factor(region),
    x = xy.593[,1],
    y = xy.593[,2],
    ) %>% # add the x and y coordinates of the genes 
    # add their log2FoldChange values to the data frame
    dplyr::mutate(log2FoldChange = rowMeans(upset.593[,5:7], na.rm = T))

  ### ------------------------- Patient 6---------------------------------------
  venn.673 <- make_vennbase(results.673$`2DN_vs_2DH`$sig_df, 
                            results.673$`2DN_vs_3D`$sig_df, 
                            results.673$`2DN_vs_Tumor`$sig_df,
                            c("2DH","3D", "Tumor"))
  # get the membership of the genes in the regions
  upset.673 <- make_upsetbase(venn.673$table, c("2DH","3D", "Tumor"))
  # get the centrum of the regions
  arranged.673 <- arrange_venn(upset.673, c("2DH","3D", "Tumor"),
                               extract_regions = T)
  # randomize the position of the genes in the regions
  set.seed(42)
  xy.673 = rbind(
    calculateCircle(x = 0.00, y = 0.2886, r = 0.1, # 3-way intersection
                    noiseFun = function(x) (x + rnorm(1,0,0.1)), # add noise
                    steps = 9,randomDist = T, randomFun = rnorm), # 2 genes
    calculateEllipse(x = -0.65, y = -0.0866, a = 0.2, b = .2, angle = -240, 
                     noiseFun = function(x) (x + rnorm(1,0,0.1)),
                     steps = 14, randomDist = T, randomFun = rnorm),
    calculateEllipse(x = 0.65, y = -0.0866, a = 0.2, b = .2, angle = -120,
                     noiseFun = function(x) (x + rnorm(1,0,0.1)),
                     steps = 217, randomDist = T, randomFun = rnorm),
    calculateEllipse(x = 0, y = 1.0392, a = 0.1, b = .2, angle = 0,
                     noiseFun = function(x) (x + rnorm(1,0,0.1)),
                     steps = 4, randomDist = T, randomFun = rnorm),
    calculateCircle(x = 0.00, y = -1.2124, r = .3,
                    noiseFun = function(x) (x + rnorm(1,0,0.2)),
                    steps = 528, randomDist = T, randomFun = rnorm),
    calculateCircle(x = 1.30, y = 1.0392, r = .3,
                    noiseFun = function(x) (x + rnorm(1,0,0.2)),
                    steps = 562,randomDist = T, randomFun = rnorm),
    calculateCircle(x = -1.30, y = 1.0392, r = .3,
                    noiseFun = function(x) (x + rnorm(1,0,0.2)),
                    steps = 59,randomDist = T, randomFun = rnorm)
  )
  # create data table for plotting
  upset.673 <- upset.673 %>%
    dplyr::mutate(region = dplyr::case_when(
      `2DH` & `3D` & `Tumor` ~ "2DH-3D-Tumor",
      `2DH` & !`3D` & `Tumor` ~ "2DH-Tumor",
      !`2DH` & `3D` & `Tumor` ~ "3D-Tumor",
      `2DH` & `3D` & !`Tumor` ~ "2DH-3D",
      !`2DH` & !`3D` & `Tumor` ~ "Tumor",
      !`2DH` & `3D` & !`Tumor` ~ "3D",
      `2DH` & !`3D` & !`Tumor` ~ "2DH"
    ), # relabel the regions
    region = as.factor(region),
    x = xy.673[,1],
    y = xy.673[,2],
    ) %>% # add the x and y coordinates of the genes 
    # add their log2FoldChange values to the data frame
    dplyr::mutate(log2FoldChange = rowMeans(upset.673[,5:7], na.rm = T))
}
# ---------------------------------------------------------------------------- #
# -                             SURFME filter                                - #
# ---------------------------------------------------------------------------- #
###  Molecular functions
con <- GO.db::GO_dbconn()
go_term <- tbl(con, "go_term")
go_mf_parent <- tbl(con, "go_mf_parents")

surfme_mf <- AnnotationDbi::select(org.Hs.eg.db, surfme$geneID, "GO","SYMBOL") %>%
  dplyr::filter(ONTOLOGY =="MF") %>%
  dplyr::mutate(EVIDENCE = factor(EVIDENCE, levels = c("EXP","HDA","IDA","IMP","IGI","TAS","IC",
                             "IEA","IBA","ISS","ND","RCA","NAS","IPI"))) %>% 
  dplyr::arrange(SYMBOL, EVIDENCE) %>% 
  dplyr::distinct(SYMBOL, .keep_all = T) %>% 
  dplyr::mutate(
    AnnotationDbi::select(GO.db, GO, c("TERM"), "GOID")
  ) %>%
  dplyr::select(!c("ONTOLOGY","GOID","EVIDENCE"))

surfme_mf <- dplyr::full_join(surfme, surfme_mf, by = c("geneID"="SYMBOL")) %>%
  dplyr::mutate(
    GOID = ifelse(is.na(GO), "GO:0003674", GO),
    TERM = ifelse(is.na(TERM), "molecular_function", TERM)) %>% 
  dplyr::select(!GO)
  

if (exists("SURFME.593") == F) {
  # 1.) Extract the significant genes from the differential expression analysis
  ### ------- Patient 593 --------
  SURFME.593 <- list()
  SURFME.593$list <- list(
    results.593$`2DN_vs_2DH`$sig_df[,c("geneID","log2FoldChange")],
    results.593$`2DN_vs_3D`$sig_df[,c("geneID","log2FoldChange")],
    results.593$`2DN_vs_Tumor`$sig_df[,c("geneID","log2FoldChange")]) 
  
  SURFME.593$upset_base <- SURFME.593$list %>% 
    merge.rec(., by="geneID", all=T, suffixes=c("", "")) %>% 
    setNames(., c("geneID","2DN_vs_2DH","2DN_vs_3D","2DN_vs_Tumor")) %>%
    # trend in regulation
    rowwise() %>%
    mutate(
      positivity = case_when(
        all(c_across(`2DN_vs_2DH`:`2DN_vs_Tumor`)[!is.na(c_across(`2DN_vs_2DH`:`2DN_vs_Tumor`))] > 0) ~ "up-reg",
        all(c_across(`2DN_vs_2DH`:`2DN_vs_Tumor`)[!is.na(c_across(`2DN_vs_2DH`:`2DN_vs_Tumor`))] < 0) ~ "down-reg",
        TRUE ~ "mixed"
      )
    ) %>% 
    dplyr::ungroup() %>% 
    merge(., surfme[,c("geneID")], by="geneID") %>% 
    dplyr::mutate(across(`2DN_vs_2DH`:`2DN_vs_Tumor`, ~ifelse(is.na(.), F, T)))
  
  ### Merge with molecular functions
  SURFME.593$upset_base <- merge(SURFME.593$upset_base, surfme_mf[,c("geneID", "GOID", "TERM")], 
                                 by="geneID", all.x = T)
  SURFME.593$upset_base <- map_to_super_family(.df = SURFME.593$upset_base,
                                               .ref = surfme_mf$GOID)
  
  SURFME.593$upset_base$TERM.parent <- factor(SURFME.593$upset_base$TERM.parent,
                                              levels = na.omit(levels(SURFME.593$upset_base$TERM.parent)[match(names(MF_colors), 
                                                                                         levels(SURFME.593$upset_base$TERM.parent))]))
  
  ### ------- Patient 673 --------
  SURFME.673 <- list()
  SURFME.673$list <- list(
    results.673$`2DN_vs_2DH`$sig_df[,c("geneID","log2FoldChange")],
    results.673$`2DN_vs_3D`$sig_df[,c("geneID","log2FoldChange")],
    results.673$`2DN_vs_Tumor`$sig_df[,c("geneID","log2FoldChange")]) 
  
  SURFME.673$upset_base <- SURFME.673$list %>% 
    merge.rec(., by="geneID", all=T, suffixes=c("", "")) %>% 
    setNames(., c("geneID","2DN_vs_2DH","2DN_vs_3D","2DN_vs_Tumor")) %>%
    # trend in regulation
    rowwise() %>%
    mutate(
      positivity = case_when(
        all(c_across(`2DN_vs_2DH`:`2DN_vs_Tumor`)[!is.na(c_across(`2DN_vs_2DH`:`2DN_vs_Tumor`))] > 0) ~ "up-reg",
        all(c_across(`2DN_vs_2DH`:`2DN_vs_Tumor`)[!is.na(c_across(`2DN_vs_2DH`:`2DN_vs_Tumor`))] < 0) ~ "down-reg",
        TRUE ~ "mixed"
      )
    ) %>% 
    dplyr::ungroup() %>% 
    merge(., surfme[,c("geneID")], by="geneID") %>% 
    dplyr::mutate(across(`2DN_vs_2DH`:`2DN_vs_Tumor`, ~ifelse(is.na(.), F, T)))
  
  ### Merge with molecular functions
  SURFME.673$upset_base <- merge(SURFME.673$upset_base, surfme_mf[,c("geneID", "GOID", "TERM")], 
                                 by="geneID", all.x = T)
  SURFME.673$upset_base <- map_to_super_family(.df = SURFME.673$upset_base,
                                               .ref = surfme_mf$GOID)
  SURFME.673$upset_base$TERM.parent <- factor(SURFME.673$upset_base$TERM.parent,
                                              levels = levels(
                                                SURFME.673$upset_base$TERM.parent)[match(names(MF_colors), 
                                                                                         levels(SURFME.673$upset_base$TERM.parent))])
  }


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



# Load proteomics results
imputed_proteomics <- readxl::read_excel(file.path("data", "WTF2imputed.xlsx"))
proteomics_IDs <- imputed_proteomics[,c("Accession","Gene Symbol")]
proteomics_IDs <- setNames(proteomics_IDs, c("UniprotID","geneID"))

proteomics_matrix <- imputed_proteomics %>% 
  dplyr::select(-c(1:4)) %>% 
  dplyr::select(contains(c("2D Normoxia","3D","Tumor"))) %>% 
  dplyr::select(contains(c("593","673"))) %>% 
  dplyr::select(contains(c("Surface"))) %>% 
  setNames(., c("593_2D","593_3D","593_Tumor","673_2D","673_3D","673_Tumor")) %>% 
  as.matrix()

rownames(proteomics_matrix) <- proteomics_IDs$UniprotID
proteomics_matrix <- proteomics_matrix[which(complete.cases(proteomics_matrix)),]

proteomics_IDs <- proteomics_IDs[which(proteomics_IDs$UniprotID %in% 
                                         surfme$UniprotID),]

protein_surfme <- proteomics_matrix[which(rownames(proteomics_matrix) %in% 
                                           proteomics_IDs$UniprotID),]
protein_surfme <- normalize.quantiles(protein_surfme, keep.names = T)
