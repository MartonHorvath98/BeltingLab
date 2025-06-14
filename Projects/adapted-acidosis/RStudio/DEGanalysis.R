# Created: 2024 10  07 ; Last Modified: 2025 01 20
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

DEG_palette <- c("NS" = '#c1c1c1',
                 "Log10P" = '#363636',
                 "Log2FoldChange" = '#767676',
                 "Signif. up-regulated" = '#841f27',
                 "Signif. down-regulated" = '#000f64')

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
    stringr::str_detect(samples, "sel") ~ "acidosis",
    stringr::str_detect(samples, "acu") ~ "acidosis",
    stringr::str_detect(samples, "hypo") ~ "hypoxia"
  ),
  treatment = factor(treatment, 
                     levels = c("control","hypoxia",
                                "acidosis")))

###  Data exploration
U87.plots <- list()
# create PCA plot
pca_base <- prcomp(t(U87.expr[,5:25]), center = TRUE, scale. = TRUE)
(U87.plots$PCA <- plot_pca(
  data = pca_base, 
  .groups = U87.meta$samplenames, 
  .labels = c(
    "control_sel"="Control (CA)","sel_pH64"="Chronic acidosis (pH 6.4)",
    "control_acu"="Control (AA)","acu_pH64"="Acute acidosis (pH 6.4)","acu_pH68"="Acute acidosis (pH 6.8)",
    "control_nox"="Control (NOX)","hypoxia"="Hypoxia"),
  .values = c(1:7)))

(U87.plots$PCA <- U87.plots$PCA +
    # define legend categories
    scale_shape_manual(name = "Groups",
                       labels = c(
                         "control_sel"="Control (CA)","sel_pH64"="Chronic acidosis (pH 6.4)",
                         "control_acu"="Control (AA)","acu_pH64"="Acute acidosis (pH 6.4)","acu_pH68"="Acute acidosis (pH 6.8)",
                         "control_nox"="Control (NOX)","hypoxia"="Hypoxia"),
                       values = c(
                         "control_sel"=15,"sel_pH64"=15,
                         "control_acu"=15,"acu_pH64"=15,"acu_pH68"=15,
                         "control_nox"=15,"hypoxia"=15)) +
  guides(color=guide_legend(ncol=3), fill=guide_legend(ncol=3)))

### Over-representation (ORA) and Gene Set Enrichment Analysis (GSEA)
U87.deg <- lapply(U87.deg, function(x){
  df = x %>% 
    dplyr::rename(PROBEID = "ID.ID", Symbol = "ID.Symbol",
                  entrezID = "ID.entrezID", log2FoldChange = "logFC",
                  pvalue = "P.Value", padj = "adj.P.Val")
  
  get_significance(.df = df)
})
# extract entrez IDs for gene set of interest and background
U87.genes <- list()
U87.genes <- lapply(U87.deg, function(x){
  get_genelist(.df = x, 
               .filter = x[["significance"]] %in% c("Signif. up-regulated", 
                                                    "Signif. down-regulated"),
               .value = "t",
               .name = "entrezID")
})
# Run the ORA analysis
U87.ORA <- list()
U87.ORA <- lapply(U87.genes,
                   function (x){
                     run_ora(.interest = x[["interest"]],
                             .background = x[["background"]],
                             .pathways = pathways)})
U87.ORA <- lapply(U87.ORA, function(x){
  return(c(x, extract_ora_results(.ora = x[["ora"]],
                                  .db = pathways)))})
# Save the results
dir.create(file.path(wd, results_dir, date, "tables"), recursive = T, showWarnings = FALSE)
sapply(names(U87.ORA), function(x){
  openxlsx::write.xlsx(U87.ORA[[x]][c(2:3)], 
                       file.path(wd, results_dir, date, "tables",
                                 paste0("U87_", x, "_all_pathway_ORA.xlsx")))
})

# Run the GSEA analysis
U87.GO <- list()
U87.GO <- lapply(U87.genes, function(x){
  run_gsea(.geneset = x[["background"]], .terms = terms)})
U87.GO <- lapply(U87.GO,
                  function(x){
                    return(c(x, 
                             extract_gsea_results(
                               .gsea = x[["gsea"]],
                               .db = terms)))
                  })
# Save the results
sapply(names(U87.GO), function(x){
  openxlsx::write.xlsx(U87.GO[[x]][c(2:3)], 
                       file.path(wd, results_dir, date, "tables",
                                 paste0("U87_", x, "_all_go&hallmark_GSEA.xlsx")))
})

U87.GSEA <- list()
U87.GSEA <- lapply(U87.genes, function(x){
  run_gsea(.geneset = x[["background"]], .terms = pathways)})
U87.GSEA <- lapply(U87.GSEA,
                    function(x){
                      return(c(x, 
                               extract_gsea_results(
                                 .gsea = x[["gsea"]],
                                 .db = pathways)))
                    })
# Save the results
sapply(names(U87.GSEA), function(x){
  openxlsx::write.xlsx(U87.GSEA[[x]][c(2:3)], 
                       file.path(wd, results_dir, date, "tables",
                                 paste0("U87_", x, "_all_pathway_GSEA.xlsx")))
})


interest_genes_U87 <- c("CSGALNACT1","CSGALNACT2","CHSY1","CHSY3","CHPF","CHPF2",
                        "DCN","BGN","DSE","DSEL","CA9","HIF1A","G0S2","C7orf68")

U87.subset <- list()

U87.subset$logexp <- U87.expr %>%
  dplyr::select(!ENTREZID) %>% 
  dplyr::filter(SYMBOL %in% interest_genes & 
                  PROBEID %in% U87.deg$`sel_pH647-control_sel`$ID.ID) 

U87.subset$exp <- U87.subset$logexp %>% 
  dplyr::mutate(across("200118400068_I":"200118400033_I", ~ 2^.))

U87.subset$df <- dplyr::inner_join(
  U87.subset$exp %>% 
    tidyr::pivot_longer(cols = !c(SYMBOL,PROBEID), values_to = "expression",
                        names_to = "sample"),
  U87.subset$logexp %>% 
    tidyr::pivot_longer(cols = !c(SYMBOL,PROBEID), values_to = "logexp",
                        names_to = "sample"),
  by = c("PROBEID","SYMBOL","sample")) %>% 
  dplyr::rowwise(.) %>% 
  dplyr::mutate(treatment = case_match(sample,
    c("200118400068_I","200118400035_I","200118400035_D") ~ "control_sel",
    c("200118400035_K","200118400035_G","200118400035_F") ~ "sel_pH64",
    c("200118400068_K","200118400068_D","200118400068_A") ~ "control_acu",
    c("200118400068_C","200118400068_G","200118400033_E") ~ "acu_pH68",
    c("200118400033_B","200118400068_B","200118400035_B") ~ "acu_pH64",
    c("200118400033_G","200118400035_L","200118400033_H") ~ "control_nox",
    c("200118400035_E","200118400035_H","200118400033_I") ~ "hypoxia"
  )) %>%
  dplyr::group_by(SYMBOL, treatment) %>% 
  dplyr::reframe(sample = sample, 
                 expression = expression,
                 mean = mean(expression),
                 se = sd(expression)/sqrt(n()),
                 logexp = logexp)


tmp <- limmaDEA(.data = U87.subset$logexp,
         .design = samples,
         .contrast = c("sel_pH647-control_sel",
                       "acu_pH68-control_acu",
                       "acu_pH64-control_acu",
                       "hypoxia-control_nox"))

tmp <- setNames(tmp, c("sel_pH64","acu_pH68",
                "acu_pH64","hypoxia"))

lapply(names(tmp), function(x){
  return(tmp[[x]] %>% 
           dplyr::mutate(
             treatment = x,
             SYMBOL = ID.Symbol) %>% 
           dplyr::select(SYMBOL,treatment, logFC))
}) %>% rbind.fill(.) %>% 
  dplyr::full_join(U87.subset$df, ., by = c("SYMBOL", "treatment")) %>% 
  dplyr::mutate(logFC = ifelse(is.na(logFC), 0, logFC)) -> U87.subset$df



dir.create(file.path(wd, "Results", "U87", date, "tables"), recursive = T, showWarnings = FALSE)
openxlsx::write.xlsx(U87.subset, 
                     file.path(wd, "Results", "U87", date, "tables",
                               "U87_linear_scale_genes-of-interest_expression.xlsx"))


# ---- 4.) PANC1 Chronic Acidosis AA vs NA (Clariom D Human Pico) Affymetrix ---- #
results_dir <- "Results/PANC1"
dir.create(file.path(wd, results_dir, date), recursive = T, showWarnings = FALSE)

PANC1.deg <- PANC1.deg %>% 
  dplyr::rename(PROBEID = "ID.ID", Symbol = "ID.Symbol",
                entrezID = "ID.entrezID", log2FoldChange = "logFC",
                pvalue = "P.Value", padj = "adj.P.Val") %>% 
  get_significance(.)

table(PANC1.deg$significance)
# Save file
dir.create(file.path(wd, results_dir, date, "tables"), recursive = T, showWarnings = FALSE)
write.xlsx(PANC1.deg, file.path(wd, results_dir, date, "tables", "PANC1_DEGs.xlsx"))

## extract entrez IDs for gene set of interest and background
PANC1.genes <- get_genelist(.df = PANC1.deg,
                           .filter = PANC1.deg[["significance"]] %in% c("Signif. up-regulated", 
                                                                      "Signif. down-regulated"),
                           .value = "log2FoldChange",
                           .name = "entrezID")

## Run the ORA analysis
PANC1.ORA <- run_ora(.interest = PANC1.genes$interest,
                    .background = PANC1.genes$background,
                    .pathways = pathways)
PANC1.ORA <- c(PANC1.ORA, extract_ora_results(.ora = PANC1.ORA$ora, .db = pathways))
# Save the results
openxlsx::write.xlsx(PANC1.ORA[c(2:3)], 
                     file.path(wd, results_dir, date, "tables", "PANC1_pathway_ORA.xlsx"))

# GSEA on the GO terms and hallmarks from MSigDB
PANC1.GO <- run_gsea(.geneset = PANC1.genes$background, 
                    .terms = terms)
PANC1.GO <- c(PANC1.GO, extract_gsea_results(.gsea = PANC1.GO$gsea, .db = terms))
# Save the results
openxlsx::write.xlsx(PANC1.GO[c(2:3)], 
                     file.path(wd, results_dir, date, "tables", "PANC1_all_go&hallmark_GSEA.xlsx"))

# GSEA on the KEGG- and REACTOME pathways from MSigDB
PANC1.GSEA <- run_gsea(.geneset = PANC1.genes$background,
                      .terms = pathways)
PANC1.GSEA <- c(PANC1.GSEA, extract_gsea_results(.gsea = PANC1.GSEA$gsea, .db = pathways))
# Save the results
openxlsx::write.xlsx(PANC1.GSEA[c(2:3)],
                     file.path(wd, results_dir, date, "tables", "PANC1_all_pathway_GSEA.xlsx"))


interest_genes_PANC1 <- c("CSGALNACT1","CSGALNACT2","CHSY1","DCN","BGN","DSE","DSEL")

PANC1.subset <- list()

PANC1.subset$logexp <- PANC1.expr %>%
  dplyr::select(!c(ENTREZID, GENENAME)) %>% 
  dplyr::filter(SYMBOL %in% interest_genes_PANC1 & 
                  PROBEID %in% PANC1.deg$PROBEID) 
                

PANC1.subset$exp <- PANC1.subset$logexp %>% 
  dplyr::mutate(across("1_NA1_(Clariom_D_Human).CEL":"6_AA3_(Clariom_D_Human).CEL", ~ 2^.))

PANC1.subset$df <- dplyr::inner_join(
  PANC1.subset$exp %>% 
    tidyr::pivot_longer(cols = !c(SYMBOL,PROBEID), values_to = "expression",
                        names_to = "sample"),
  PANC1.subset$logexp %>% 
    tidyr::pivot_longer(cols = !c(SYMBOL,PROBEID), values_to = "logexp",
                        names_to = "sample"),
  by = c("PROBEID","SYMBOL","sample")) %>% 
  dplyr::rowwise(.) %>% 
  dplyr::mutate(treatment = case_match(sample,
                                       c("1_NA1_(Clariom_D_Human).CEL",
                                         "2_NA2_(Clariom_D_Human).CEL",
                                         "3_NA3_(Clariom_D_Human).CEL") ~ "PANC1_NA",
                                       c("4_AA1_(Clariom_D_Human).CEL",
                                         "5_AA2_(Clariom_D_Human).CEL",
                                         "6_AA3_(Clariom_D_Human).CEL") ~ "PANC1_AA"
  )) %>%
  dplyr::group_by(SYMBOL, treatment) %>% 
  dplyr::reframe(sample = sample, 
                 expression = expression,
                 mean = mean(expression),
                 se = sd(expression)/sqrt(n()),
                 logexp = logexp)


tmp <- limmaDEA(.data = PANC1.subset$logexp,
                .design = c("PANC1_NA", "PANC1_NA","PANC1_NA",
                            "PANC1_AA", "PANC1_AA", "PANC1_AA"),
                .contrast = c("PANC1_AA - PANC1_NA"))

tmp[[1]] %>% 
  dplyr::mutate(
    SYMBOL = ID.Symbol) %>% 
  dplyr::select(SYMBOL, logFC) %>% 
  rbind.fill(.) %>% 
  dplyr::full_join(PANC1.subset$df, ., by = c("SYMBOL")) %>% 
  dplyr::mutate(logFC = ifelse(is.na(logFC), 0, logFC)) -> PANC1.subset$df



dir.create(file.path(wd, "Results", "PANC1", date, "tables"), recursive = T, showWarnings = FALSE)
openxlsx::write.xlsx(PANC1.subset, 
                     file.path(wd, "Results", "PANC1", date, "tables",
                               "PANC1_linear_scale_genes-of-interest_expression.xlsx"))

# 
# tmp <- merge.rec(list(CCLD.df[,c(1,2,4)] %>% dplyr::filter(Symbol %in% interest_genes$gene_symbol),
#                HGCC.deg$U3017[,c(2,4,8,10)] %>% dplyr::filter(Symbol %in% interest_genes$gene_symbol),
#                HGCC.deg$U3047[,c(2,4,8,10)] %>% dplyr::filter(Symbol %in% interest_genes$gene_symbol),
#                HGCC.deg$U3054[,c(2,4,8,10)] %>% dplyr::filter(Symbol %in% interest_genes$gene_symbol)),
#                #U87.deg$`sel_pH647-control_sel`[,c(2,10)] %>% dplyr::filter(Symbol %in% interest_genes$Gene.name),
#                #PANC1.deg[,c(2,10)] %>% dplyr::filter(Symbol %in% interest_genes$Gene.name)),
#           by = "Symbol", all = T, suffixes = c("",""))
# tmp <- setNames(tmp, c("Symbol","log2FC.CCLD","regulation.CCLD",
#                        "log2FC.U3017","padj.U3017","regulation.U3017",
#                        "log2FC.U3047","padj.U3047","regulation.U3047",
#                        "log2FC.U3054","padj.U3054","regulation.U3054"))
# interest_genes <- merge(interest_genes, tmp, by.x = "gene_symbol", by.y = "Symbol")                 
# write.xlsx(interest_genes, "Results/selected-5-category-genes.xlsx")
