# Created: 2024 10  07 ; Last Modified: 2025 01 16
# MH
################################################################################
#                     Set up the environment                                   #
################################################################################
# Set working directory
wd <- getwd()

path <- c("data/processed/")
files <- list.files(file.path(wd, path),
                    pattern = ".xlsx$", full.names = TRUE)
date <- Sys.Date()
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
dir.create(file.path(wd, results_dir, date), recursive = T, showWarnings = FALSE)

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
dir.create(file.path(wd, results_dir, date, "tables"), recursive = T, showWarnings = FALSE)
write.xlsx(merge(CCLD.df, CCLD.expr[,c(2,3,5)], by.x = "Symbol", by.y = "SYMBOL"),
           file.path(wd, results_dir, date, "tables", "LDvsnoLD_DEGs.xlsx"))

## extract entrez IDs for gene set of interest and background
CCLD.genes <- get_genelist(.df = CCLD.df,
                           .filter = CCLD.df[["significance"]] %in% c("Signif. up-regulated", 
                                                                "Signif. down-regulated"),
                           .value = "log2FoldChange",
                           .name = "entrezID")
## Run the ORA analysis
CCLD.ORA <- list()
CCLD.ORA <- run_ora(.interest = CCLD.genes$interest,
                    .background = CCLD.genes$background,
                    .pathways = pathways)

CCLD.ORA <- c(CCLD.ORA, extract_ora_results(.ora = CCLD.ORA$ora, .db = pathways))
# Save the results
openxlsx::write.xlsx(CCLD.ORA[c(2:3)], 
                     file.path(wd, results_dir,  date, "tables", "LDvsnoLD_pathway_ORA.xlsx"))

# GSEA on the GO terms from MSigDB
CCLD.GO <- list()
CCLD.GO <- run_gsea(.geneset = CCLD.genes$background, 
                    .terms = terms)
CCLD.GO <- c(CCLD.GO, extract_gsea_results(.gsea = CCLD.GO$gsea, .db = terms))

# Save the results
openxlsx::write.xlsx(CCLD.GO[c(2:3)], 
                     file.path(wd, results_dir,  date, "tables", "LDvsnoLD_all_go&hallmark_GSEA.xlsx"))

# GSEA on the KEGG- and REACTOME pathways from MSigDB
CCLD.GSEA <- list()
CCLD.GSEA <- run_gsea(.geneset = CCLD.genes$background,
                      .terms = pathways)
CCLD.GSEA <- c(CCLD.GSEA, extract_gsea_results(.gsea = CCLD.GSEA$gsea, .db = pathways))

# Save the results
openxlsx::write.xlsx(CCLD.GSEA[c(2:3)],
                     file.path(wd, results_dir,  date, "tables", "LDvsnoLD_all_pathway_GSEA.xlsx"))

# ---- 2.) Hugo's Primary cells 2D vs 3D (Clariom D Human Pico) Affymetrix ---- #
results_dir <- "Results/HGCC"
dir.create(file.path(wd, results_dir, date), recursive = T, showWarnings = FALSE)

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
HGCC.deg <- lapply(HGCC.deg, function(x){
  df = x %>% 
    dplyr::rename(PROBEID = "ID.ID", Symbol = "ID.Symbol", 
                  entrezID = "ID.Entrez", log2FoldChange = "logFC",
                  pvalue = "P.Value", padj = "adj.P.Val")
  
  get_significance(.df = df)
})

HGCC.venn <- list()
HGCC.venn <- get_regions(.list = HGCC.deg[1:3], .names = c("U3017","U3047","U3054"))

(HGCC.plots$VENN <- plot_venn(.data = HGCC.venn$df,
                           .sets = c("U3017","U3047","U3054"),
                           .labels = c("U3054",
                                       "U3047",
                                       "U3017")))

# Create a Vulcano plot for the shared DEGs
(p1 <- plot_vulcan(HGCC.deg$U3017) + ggtitle("U3017"))
(p2 <- plot_vulcan(HGCC.deg$U3047) + ggtitle("U3047"))
(p3 <- plot_vulcan(HGCC.deg$U3054) + ggtitle("U3054"))

(HGCC.plots$vulcano <- cowplot::plot_grid(p1, p2, p3, nrow = 1))

dir.create(file.path(wd, results_dir, date, "plots"), recursive = T, showWarnings = FALSE)
# Save the PCA plot
ggsave(file.path(wd, results_dir, date, "plots", "HGCC_PCA_plot.png"), bg = "white",
       HGCC.plots$PCA, width = 8, height = 6, dpi = 300)
# Save the VENN plot
ggsave(file.path(wd, results_dir, date, "plots", "HGCC_VENN_plot.png"), bg = "white",
       HGCC.plots$VENN, width = 12, height = 8, dpi = 300)
# Save the Vulcano plot
ggsave(file.path(wd, results_dir, date, "plots", "HGCC_Vulcano_plot.png"), bg = "white",
       HGCC.plots$vulcano, width = 18, height = 6, dpi = 300)

rm(list = c("pca_base","p1","p2","p3","samples"))

### Over-representation (ORA) and Gene Set Enrichment Analysis (GSEA)
# extract entrez IDs for gene set of interest and background
HGCC.genes <- list()
HGCC.genes <- lapply(HGCC.deg, function(x){
  get_genelist(.df = x, 
               .filter = x[["significance"]] %in% c("Signif. up-regulated", 
                                                    "Signif. down-regulated"),
               .value = "t",
               .name = "entrezID")
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
dir.create(file.path(wd, results_dir, date, "tables"), recursive = T, showWarnings = FALSE)
sapply(names(HGCC.ORA), function(x){
  openxlsx::write.xlsx(HGCC.ORA[[x]][c(2:3)], 
                       file.path(wd, results_dir, date, "tables",
                                 paste0("HGCC_", x, "_all_pathway_ORA.xlsx")))
})

# Run the GSEA analysis
HGCC.GO <- list()
HGCC.GO <- lapply(HGCC.genes, function(x){
  run_gsea(.geneset = x[["background"]], .terms = terms)})
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
                       file.path(wd, results_dir, date, "tables",
                                 paste0("HGCC_", x, "_all_go&hallmark_GSEA.xlsx")))
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
                       file.path(wd, results_dir, date, "tables",
                                 paste0("HGCC_", x, "_all_pathway_GSEA.xlsx")))
})

# ---- 3.) U87 Chronic Acidosis AA vs NA & HOX vs NOX (Illumina BeadChip) ---- #
results_dir <- "Results/U87"
dir.create(file.path(wd, results_dir, date), recursive = T, showWarnings = FALSE)

# Prepare metadata table based on the sample names
samples <- factor(c("control_sel", "control_sel", "control_sel",
                    "sel_pH64", "sel_pH64", "sel_pH64",
                    "control_acu", "control_acu", "control_acu",
                    "acu_pH68", "acu_pH68", "acu_pH68",
                    "acu_pH64", "acu_pH64", "acu_pH64",
                    "control_nox", "control_nox", "control_nox",
                    "hypoxia", "hypoxia", "hypoxia"),
                  levels = c("control_sel", "sel_pH64", "control_acu",
                             "acu_pH68", "acu_pH64", "control_nox", "hypoxia"))
U87.meta <- data.frame("samplenames" = samples) %>%
  # extract information from the sample names
  dplyr::mutate(group = dplyr::case_when( # extract the cell IDs
    stringr::str_detect(samples, "sel") ~ "selection",
    stringr::str_detect(samples, "acu") ~ "acute",
    stringr::str_detect(samples, "ox") ~ "oxygen"
    ),
  group = factor(group, levels = c("oxygen","acute","selection"))) %>%
  dplyr::mutate(treatment = dplyr::case_when( # extract the cell IDs
    stringr::str_detect(samples, "control") ~ "control",
    stringr::str_detect(samples, "sel") ~ "selective-acidosis",
    stringr::str_detect(samples, "acu") ~ "acute-acidosis",
    stringr::str_detect(samples, "hypo") ~ "hypoxia"
  ),
  treatment = factor(treatment, 
                     levels = c("control","hypoxia",
                                "acute-acidosis","selective-acidosis")))

###  Data exploration
U87.plots <- list()
# create PCA plot
pca_base <- prcomp(t(U87.expr[,2:22]), center = TRUE, scale. = TRUE)
(U87.plots$PCA <- plot_pca(
  data = pca_base, 
  .groups = U87.meta$samplenames, 
  .labels = c(
    "control_sel"="Control (SA)","sel_pH64"="Selective acidosis (pH 6.4)",
    "control_acu"="Control (AA)","acu_pH64"="Acute acidosis (pH 6.4)","acu_pH68"="Acute acidosis (pH 6.8)",
    "control_nox"="Control (NOX)","hypoxia"="Hypoxia"),
  .values = c(
    "control_sel"="skyblue","sel_pH64"="salmon",
    "control_acu"="skyblue","acu_pH64"="magenta","acu_pH68"="purple",
    "control_nox"="skyblue","hypoxia"="darkred")))

(U87.plots$PCA <- U87.plots$PCA +
    # define legend categories
    scale_shape_manual(name = "Groups",
                       labels = c(
                         "control_sel"="Control (SA)","sel_pH64"="Selective acidosis (pH 6.4)",
                         "control_acu"="Control (AA)","acu_pH64"="Acute acidosis (pH 6.4)","acu_pH68"="Acute acidosis (pH 6.8)",
                         "control_nox"="Control (NOX)","hypoxia"="Hypoxia"),
                       values = c(
                         "control_sel"=15,"sel_pH64"=15,
                         "control_acu"=16,"acu_pH64"=16,"acu_pH68"=16,
                         "control_nox"=17,"hypoxia"=17)) +
  guides(color=guide_legend(ncol=3), fill=guide_legend(ncol=3)))


interest_genes <- c("CSGALNACT1","CSGALNACT2","CHSY1","CHSY3","CHPF","CHPF2",
                    "DCN","BGN","DSE","DSEL","CA9","HIF1A","G0S2","C7orf68")

U87.subset <- as.data.frame(U87.dat) %>% 
  dplyr::select(c("SYMBOL" | contains(samples))) %>% 
  tibble::rownames_to_column("PROBEID") %>% 
  dplyr::filter(PROBEID %in% c(U87.deg$`sel_pH647-control_sel`$ID.ID,
                                      U87.deg$`acu_pH68-control_acu`$ID.ID,
                                      U87.deg$`acu_pH64-control_acu`$ID.ID,
                                      U87.deg$`hypoxia-control_nox`$ID.ID)) %>% 
  dplyr::filter(SYMBOL %in% interest_genes)

U87.subset <- U87.subset %>% 
  tidyr::pivot_longer(cols = !c(SYMBOL,PROBEID), values_to = "expression",
                      names_to = "sample") %>% 
  dplyr::rowwise(.) %>% 
  dplyr::mutate(treatment = case_match(sample,
    c("200118400068_I","200118400035_I","200118400035_D") ~ "control_sel",
    c("200118400035_K","200118400035_G","200118400035_F") ~ "sel_pH64",
    c("200118400068_K","200118400068_D","200118400068_A") ~ "control_acu",
    c("200118400068_C","200118400068_G","200118400033_E") ~ "acu_pH64",
    c("200118400033_B","200118400068_B","200118400035_B") ~ "acu_pH68",
    c("200118400033_G","200118400035_L","200118400033_H") ~ "control_nox",
    c("200118400035_E","200118400035_H","200118400033_I") ~ "hypoxia"
  )) %>%
  dplyr::group_by(SYMBOL, treatment) %>% 
  dplyr::summarize(sample = sample, 
                   expression = expression,
                   mean = mean(expression),
                   se = sd(expression)/sqrt(n()))

dir.create(file.path(wd, results_dir, date, "tables"), recursive = T, showWarnings = FALSE)
openxlsx::write.xlsx(U87.subset, 
                     file.path(wd, results_dir, date, "tables",
                               "U87_linear_scale_genes-of-interest_expression.xlsx"))
