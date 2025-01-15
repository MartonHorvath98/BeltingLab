# Created: 2024 10  07 ; Last Modified: 2024 10 21
# MH

################################################################################
#                     Set up the environment                                   #
################################################################################
# Set working directory
wd <- getwd()

path <- c("data/processed/")
files <- list.files(file.path(wd, path),
                    pattern = ".xlsx$", full.names = TRUE)
# Load packages
source(file = file.path(wd, "packages.R"))
# Source data processing functions
source(file = file.path(wd, "functions.R"))

# load databases
msigdbr_df <- msigdbr(species = "Homo sapiens") # save local database

pathways <- msigdbr_df %>% 
  dplyr::filter(
    gs_cat == "C2", # only canonical representations (compiled by experts)
    gs_subcat %in% c("CP:KEGG", "CP:REACTOME") # KEGG and Reactome pathways
  )

terms <- msigdbr_df %>% 
  dplyr::filter(
    gs_cat == "H" | # Hallmark gene sets
    gs_cat == "C5" & # Ontology gene sets
    gs_subcat == "GO:BP" # GO terms, biological processes
  )
terms <- terms %>% 
  dplyr::rowwise(.) %>% 
  dplyr::mutate(gs_subcat = ifelse(gs_cat == "H", "HALLMARK", gs_subcat),
                gs_exact_source = ifelse(gs_cat == "H", gs_id, gs_exact_source)) %>% 
  dplyr::ungroup(.)

################################################################################
#             Differential gene expression analysis                            #
################################################################################

# ---------------------------------------------------------------------------- #
# -    1.) Svenja's CC +LD and CC no LD (Clariom D Human Pico) Affymetrix    - #
# ---------------------------------------------------------------------------- #

#Results directory
results_dir <- "Results/CCLD"
dir.create(file.path(wd, results_dir), recursive = T, showWarnings = FALSE)

# Data exploration 
CCLD.df <- CCLD.expr[,c(5,7)] %>% 
  setNames(., c("Symbol", "log2FoldChange")) %>%
  # Add a columns...
  dplyr::mutate(
    # ... for ENTREZ gene identifiers, for downstream analyses
    entrezID = mapIds(org.Hs.eg.db, Symbol, keytype = "SYMBOL", 
                      column = "ENTREZID"),
    # ... for significance levels using the thresholds:
    #                            p-adj < 0.05, abs(log2FC) > 0.5
    significance = dplyr::case_when(
      log2FoldChange < (-1)*0.5 ~ 'Signif. down-regulated',
      log2FoldChange > 0.5 ~ 'Signif. up-regulated',
      T ~ 'NS')) %>% 
  dplyr::filter(complete.cases(.))
# Save file
write.xlsx(merge(CCLD.df, CCLD.expr[,c(2,3,5)], by.x = "Symbol", by.y = "SYMBOL"),
           file.path(wd, results_dir, "LDvsnoLD_DEGs.xlsx"))

## extract entrez IDs for gene set of interest and background
CCLD.genes <- list(
  background = CCLD.df %>%
    dplyr::arrange(desc(log2FoldChange)) %>%
    dplyr::distinct(entrezID, .keep_all = T) %>%
    dplyr::pull("log2FoldChange", name ="entrezID"),
  interest = CCLD.df %>%
    dplyr::filter(significance != "NS") %>%
    dplyr::arrange(desc(log2FoldChange)) %>%
    dplyr::pull("log2FoldChange", name ="entrezID")
)

## Run the ORA analysis
CCLD.ORA <- list()
CCLD.ORA <- run_ora(.interest = CCLD.genes$interest,
                    .background = CCLD.genes$background,
                    .pathways = pathways)

CCLD.ORA <- c(CCLD.ORA, extract_ora_results(.ora = CCLD.ORA$ora, .db = pathways))
# Save the results
openxlsx::write.xlsx(CCLD.ORA[c(2:3)], 
                     file.path(wd, results_dir, "LDvsnoLD_pathway_ORA.xlsx"))

# GSEA on the GO terms from MSigDB
CCLD.GO <- list()
CCLD.GO <- run_gsea(.geneset = CCLD.genes$background, 
                    .terms = terms)

CCLD.GO <- c(CCLD.GO, extract_gsea_results(.gsea = CCLD.GO$gsea, .db = terms))
# Save the results
openxlsx::write.xlsx(CCLD.GO[c(2:3)], 
                     file.path(wd, results_dir, "LDvsnoLD_all_go&hallmark_GSEA.xlsx"))

# GSEA on the KEGG- and REACTOME pathways from MSigDB
CCLD.GSEA <- list()
CCLD.GSEA <- run_gsea(.geneset = CCLD.genes$background,
                      .terms = pathways)
CCLD.GSEA <- c(CCLD.GSEA, extract_gsea_results(.gsea = CCLD.GSEA$gsea, .db = pathways))
openxlsx::write.xlsx(CCLD.GSEA[c(2:3)], file.path(wd, results_dir, "LDvsnoLD_all_pathway_GSEA.xlsx"))

# ---- 2.) Hugo's Primary cells 2D vs 3D (Clariom D Human Pico) Affymetrix ---- #
results_dir <- "Results/HGCC"
dir.create(file.path(wd, results_dir), recursive = T, showWarnings = FALSE)

# Prepare metadata table based on the sample names
samples <- factor(c("U3017_2D", "U3017_2D","U3017_2D",
                    "U3017_3D", "U3017_3D", "U3017_3D",
                    "U3047_2D", "U3047_2D","U3047_2D",
                    "U3047_3D","U3047_3D","U3047_3D",
                    "U3054_2D", "U3054_2D", "U3054_2D",
                    "U3054_3D", "U3054_3D","U3054_3D"))
HGCC.meta <- data.frame("samplenames" = samples) %>%
  # extract information from the sample names
  dplyr::mutate(cell_line = dplyr::case_when( # extract the cell IDs
    stringr::str_detect(samples, "U3017") ~ "U3017",
    stringr::str_detect(samples, "U3047") ~ "U3047",
    stringr::str_detect(samples, "U3054") ~ "U3054"
  ),
  cell_line = factor(cell_line, levels = c("U3017","U3047","U3054"))) %>%
  dplyr::mutate(dimension = dplyr::case_when( # extract the dimension
    stringr::str_detect(samples, "2D") ~ "2D",
    stringr::str_detect(samples, "3D") ~ "3D"
  ),
  dimension = factor(dimension, levels = c("2D","3D")))

###  Data exploration
HGCC.plots <- list()
# create PCA plot
pca_base <- prcomp(t(HGCC.expr[,2:19]), center = TRUE, scale. = TRUE)
(HGCC.plots$PCA <- plot_pca(data = pca_base, .groups = HGCC.meta$samplenames, 
                           .labels = c("U3017 (2D)", "U3017 (3D)", "U3047 (2D)", 
                                      "U3047 (3D)", "U3054 (2D)", "U3054 (3D)"),
                           .values = c("U3017_2D" = "steelblue", "U3017_3D" = "darkred",
                                      "U3047_2D" = "steelblue", "U3047_3D" = "darkred",
                                      "U3054_2D" = "steelblue", "U3054_3D" = "darkred"),
                           plot.ellipse = T, .ellipse = HGCC.meta$dimension))
  
oligo::boxplot(HGCC.expr[,2:19], col = HGCC.meta$dimension, 
               xlab = "Samples", ylab = "Expression", 
               main = "HGCC 2D vs 3D expression")
(HGCC.plots$PCA <- HGCC.plots$PCA +
  # define legend categories
  scale_fill_manual(name = "Groups",
                    values = c("2D" = "steelblue", "3D" = "darkred")) +
  scale_shape_manual(name = "Groups",
                     labels = c("U3017 (2D)", "U3017 (3D)", "U3047 (2D)", 
                                "U3047 (3D)", "U3054 (2D)", "U3054 (3D)"),
                     values = c("U3017_2D" = 15, "U3017_3D" = 15,
                                "U3047_2D" = 16, "U3047_3D" = 16,
                                "U3054_2D" = 17, "U3054_3D" = 17)) +
  guides(color=guide_legend(ncol=3), fill=guide_legend(ncol=3))) 

# Venn diagram from the shared DEGs
HGCC.df <- list(
  "U3017" = HGCC.deg[,c(3,5:7)],
  "U3047" = HGCC.deg[,c(3,8:10)],
  "U3054" = HGCC.deg[,c(3,11:13)],
  "global" = HGCC.deg[,c(3,14:16)]
)
HGCC.df <- lapply(HGCC.df, function(x){
  x %>% 
    dplyr::mutate(across(!contains("SYMBOL"), ~as.numeric(.)))
})
HGCC.df <- get_significance(.list = HGCC.df)

HGCC.venn <- list()
HGCC.venn <- get_regions(.list = HGCC.df[1:3], .names = c("U3017","U3047","U3054"))

(HGCC.plots$VENN <- plot_venn(.data = HGCC.venn$df,
                           .sets = c("U3017","U3047","U3054"),
                           .labels = c("U3054",
                                       "U3047",
                                       "U3017")))

# Create a Vulcano plot for the shared DEGs
(HGCC.plots$vulcano <- plot_vulcan(HGCC.df$global))

# Save the PCA plot
ggsave(file.path(wd, results_dir, "HGCC_PCA_plot.png"), 
       HGCC.plots$PCA, width = 8, height = 6, dpi = 300)
# Save the VENN plot
ggsave(file.path(wd, results_dir, "HGCC_VENN_plot.png"),
       HGCC.plots$VENN, width = 12, height = 8, dpi = 300)
# Save the Vulcano plot
ggsave(file.path(wd, results_dir, "HGCC_Vulcano_plot.png"),
       HGCC.plots$Vulcano, width = 8, height = 6, dpi = 300)

rm(list = c("pca_base","HGCC.meta","samples"))

### Over-representation (ORA) and Gene Set Enrichment Analysis (GSEA)
# extract entrez IDs for gene set of interest and background
HGCC.genes <- list()
HGCC.genes <- lapply(HGCC.df, function(x){
  get_genelist(.df = x, 
               .filter = x[["significance"]] %in% c("Signif. up-regulated", 
                                                    "Signif. down-regulated"))
})
# Run the ORA analysis
HGCC.ORA <- list()
HGCC.ORA <- lapply(HGCC.genes,
                   function (x){
                     run_ora(.interest = x[["interest"]],
                             .background = x[["background"]],
                             .pathways = pathways)})
HGCC.ORA <- lapply(HGCC.ORA, function(x){
  return(c(x, extract_ora_results(.ora = x[["ora"]],
                                  .db = pathways)))})
# Save the results
sapply(names(HGCC.ORA), function(x){
  openxlsx::write.xlsx(HGCC.ORA[[x]][c(2:3)], 
                       file.path(wd, results_dir, paste0("HGCC_", x, "_all_pathway_ORA.xlsx")))
})

# Run the GSEA analysis
HGCC.GO <- list()
HGCC.GO <- lapply(HGCC.genes, function(x){
  run_gsea(.geneset = x[["interest"]], .terms = terms)})
HGCC.GO <- lapply(HGCC.GO,
                  function(x){
                    return(c(x, 
                             extract_gsea_results(
                               .gsea = x[["gsea"]],
                               .db = terms)))
                  })
# Save the results
sapply(names(HGCC.GO), function(x){
  openxlsx::write.xlsx(HGCC.GO[[x]][c(2:3)], 
                       file.path(wd, results_dir, paste0("HGCC_", x, "_all_go&hallmark_GSEA.xlsx")))
})

HGCC.GSEA <- list()
HGCC.GSEA <- lapply(HGCC.genes, function(x){
  run_gsea(.geneset = x[["background"]], .terms = pathways)})
HGCC.GSEA <- lapply(HGCC.GSEA,
                    function(x){
                      return(c(x, 
                               extract_gsea_results(
                                 .gsea = x[["gsea"]],
                                 .db = pathways)))
                    })
# Save the results
sapply(names(HGCC.GSEA), function(x){
  openxlsx::write.xlsx(HGCC.GSEA[[x]][c(2:3)], 
                       file.path(wd, results_dir, paste0("HGCC_", x, "_all_pathway_GSEA.xlsx")))
})

interest_genes <- c("CSGALNACT1","CSGALNACT2","CHSY1","CHSY3","CHPF","CHPF2",
                    "DCN","BGN","DSE","DSEL","CA9","HIF1A","GOS2","HIG2")

HGCC.subset <- HGCC.deg %>% 
  dplyr::select("SYMBOL" | contains(".CEL")) %>% 
  dplyr::filter(SYMBOL %in% interest_genes)

HGCC.subset <- HGCC.subset %>% 
  tidyr::pivot_longer(cols = -SYMBOL, values_to = "expression",
                      names_to = "sample") %>% 
  dplyr::rowwise(.) %>% 
  dplyr::mutate(sample = gsub(x = sample,
                              pattern = "_parental_[123]_\\(Clariom_D_Human\\).CEL",
                              replacement = "")) %>% 
  tidyr::separate_wider_delim(col = sample, delim = "_", 
                              names = c("cell_line", "dimension")) %>%
  dplyr::group_by(SYMBOL, dimension) %>% 
  dplyr::summarize(sample = paste(cell_line, dimension,sep = "_"),
                   expression = expression,
                   mean = mean(expression),
                   se = sd(expression)/sqrt(n()))

openxlsx::write.xlsx(HGCC.subset, 
                     file.path(wd, results_dir, "HGCC_genes-of-interest_expression.xlsx"))
