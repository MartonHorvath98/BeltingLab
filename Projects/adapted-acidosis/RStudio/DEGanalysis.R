################################################################################
# Created: 2024 10 07 ; Last Modified: 2025 06 24 ; MH                         #
################################################################################
# ---------------------- Set up the environment -----------------------------  #
# Set working directory
wd <- getwd()
# Load packages
if (file.exists(file.path(wd, "packages.R"))) {
  source(file = file.path(wd, "packages.R"))
} else {
  stop("Required file 'packages.R' not found in the working directory.")
}
# Source data processing functions
if (file.exists(file.path(wd, "functions.R"))) {
  source(file = file.path(wd, "functions.R"))
} else {
  stop("Required file 'functions.R' not found in the working directory.")
}

# load databases
msigdbr_df <- msigdbr(species = "Homo sapiens") # save local database

# Filter KEGG and Reactome pathways
pathways <- msigdbr_df %>% 
  dplyr::filter(
    gs_cat == "C2", # only canonical representations (compiled by experts)
    gs_subcat %in% c("CP:KEGG", "CP:REACTOME") # KEGG and Reactome pathways
  )
# Filter Hallmakr and GO:BP gene sets
terms <- msigdbr_df %>% 
  dplyr::filter(
    gs_cat == "H" | # Hallmark gene sets
    gs_cat == "C5" & # Ontology gene sets
    gs_subcat == "GO:BP" # GO terms, biological processes
  )
# standardize gene set labels
terms <- terms %>% 
  dplyr::rowwise(.) %>% 
  dplyr::mutate(gs_subcat = ifelse(gs_cat == "H", "HALLMARK", gs_subcat),
                gs_exact_source = ifelse(gs_cat == "H", gs_id, gs_exact_source)) %>% 
  dplyr::ungroup(.)
# Save gene sets
save(msigdbr_df, pathways, terms, file = "./RData/MSigDB_gene_sets.RData")

# add colour palette
DEG_palette <- c("NS" = '#c1c1c1',
                 "Log10P" = '#363636',
                 "Log2FoldChange" = '#767676',
                 "Signif. up-regulated" = '#841f27',
                 "Signif. down-regulated" = '#000f64')

################################################################################
#             Differential gene expression analysis                            #
################################################################################
# Get current date for results directory
date <- format(Sys.Date(), "%Y-%m-%d")
# Create results directory
degs_dir <- "Results/DEG"
enrich_dir <- "Results/Enrichment"
# Create tables and plots directories
dir.create(file.path(wd, degs_dir, "tables"), recursive = T, showWarnings = FALSE)
dir.create(file.path(wd, degs_dir, "plots"), recursive = T, showWarnings = FALSE)
dir.create(file.path(wd, enrich_dir, "tables"), recursive = T, showWarnings = FALSE)
dir.create(file.path(wd, enrich_dir, "plots"), recursive = T, showWarnings = FALSE)

# ---------------------------------------------------------------------------- #
# -    1.) Svenja's CC +LD and CC no LD (Clariom D Human Pico) Affymetrix    - #
# ---------------------------------------------------------------------------- #
# Load the processed data
load(file.path(wd, "RData", "CCLD_processedData.RData"))
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

# Save the results
save(CCLD.expr, CCLD.df,
     file = file.path(wd, "RData", "CCLD_processedData.RData"))
# Save to file
write.xlsx(merge(CCLD.df, CCLD.expr[,c(2,3,5)], by.x = "Symbol", by.y = "SYMBOL"),
           file.path(wd, degs_dir, "tables", paste(date, "CCLD-DEGs.xlsx", sep = "-")))

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
                     file.path(wd, enrich_dir, "tables", paste(date, "CCLD-ORA-pathways.xlsx", sep = "-")))

# GSEA on the GO terms from MSigDB
CCLD.GO <- list()
CCLD.GO <- run_gsea(.geneset = CCLD.genes$background, 
                    .terms = terms)
CCLD.GO <- c(CCLD.GO, extract_gsea_results(.gsea = CCLD.GO$gsea, .db = terms))

# Save the results
openxlsx::write.xlsx(CCLD.GO[c(2:3)], 
                     file.path(wd, enrich_dir, "tables", paste(date, "CCLD-GSEA-go&hallmark.xlsx", sep = "-")))

# GSEA on the KEGG- and REACTOME pathways from MSigDB
CCLD.GSEA <- list()
CCLD.GSEA <- run_gsea(.geneset = CCLD.genes$background,
                      .terms = pathways)
CCLD.GSEA <- c(CCLD.GSEA, extract_gsea_results(.gsea = CCLD.GSEA$gsea, .db = pathways))

# Save the results
openxlsx::write.xlsx(CCLD.GSEA[c(2:3)],
                     file.path(wd, enrich_dir, "tables", paste(date, "CCLD-GSEA-pathways.xlsx", sep="-")))

# Visualize selected enrichment scores 
palette = c("EXTRACELLULAR_MATRIX_ORGANIZATION" = "#380186",
            "GLYCOSAMINOGLYCAN_METABOLISM" = "#8D00FA",
            "PROTEOGLYCAN_METABOLIC_PROCESS" = "#E70041",
            "CHONDROITIN_SULFATE_DERMATAN_SULFATE_METABOLISM" = "#FF08FF",
            "FAT_CELL_DIFFERENTIATION" = "#000000")
# Identify indices of selected pathways/GO terms in GSEA and GO results
pathway.ranks = which(CCLD.GSEA$df$Name %in% names(palette))
go.ranks = which(CCLD.GO$df$Name %in% names(palette))

# Extract running enrichment scores for selected gene sets
CCLD.enrichplot <- list()
CCLD.enrichplot$runningScores <- do.call(rbind, c(lapply(pathway.ranks, gsInfo, object = CCLD.GSEA$gsea),
                                                  lapply(go.ranks, gsInfo, object = CCLD.GO$gsea)))
# Clean up pathway/GO term names for plotting
CCLD.enrichplot$runningScores$Description <- factor(
  gsub(x = CCLD.enrichplot$runningScores$Description, pattern = "REACTOME_|GOBP_", replacement =  ""),
  levels = names(palette))

# Plot running enrichment scores for selected gene sets
p1 <- plotRunningScore(.df = CCLD.enrichplot$runningScores, 
                       .x = "x", .y = "runningScore", 
                       .color = "Description", .palette = palette)

# Plot gene ranks for selected gene sets, faceted by pathway/GO term
p2 <- plotGeneRank(.df = CCLD.enrichplot$runningScores, 
                   .x = "x", .facet = "Description~.",
                   .color = "Description", .palette = palette)

# Combine running score and gene rank plots into a single figure
(CCLD.enrichplot$plot <- cowplot::plot_grid(
  p1, p2 + theme(strip.text.y.left = element_blank()),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r",
  margins = c(0.5, 0.5, 0.5, 0.5)))

# Save plot as svg for publication
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "CCLD-enrichplot.svg", sep = "-")), 
       device = "svg", plot = CCLD.enrichplot$plot, 
       bg = "white", width = 6, height = 4.45, units = "in")
# Save plot as png for presentation
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "CCLD-enrichplot.png", sep = "-")), 
       device = "png", plot = CCLD.enrichplot$plot, 
       bg = "white", dpi = 300, width = 6, height = 4.45, units = "in")

# Prepare summary table of enrichment results for selected gene sets
CCLD.enrichplot$table <- rbind(CCLD.GSEA$df[pathway.ranks,],
                               CCLD.GO$df[go.ranks,]) %>% 
  dplyr::arrange(order(match(Name, names(palette)))) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("Name") %>%
  dplyr::mutate(
    NES = round(NES, 4),
    FDR = round(p.adjust, 4)) %>%
  dplyr::mutate(
    signif.level = case_when(
      # ***P< 0.01, **P< 0.05, and *P< 0.1
      FDR < 0.01 ~ paste("***"),
      FDR < 0.05 ~ paste("**"),
      FDR < 0.1 ~ paste("*"),
      TRUE ~ as.character(FDR),
    )) %>% 
  dplyr::select(c(5,13,14))

# Save table
openxlsx::write.xlsx(CCLD.enrichplot$table, 
                     file.path(wd, enrich_dir, "tables", paste(date, "CCLD-enrichplot-table.xlsx", sep = "-")))

# Save enrichment results
save(CCLD.GO, CCLD.GSEA,
     file = file.path(wd, "RData", "CCLD_enrichment_results.RData"))
# Clean up the environment
rm(list = c(ls(pattern = "CCLD.*"), ls(pattern = "*.ranks"), "p1", "p2"))

# ---------------------------------------------------------------------------- #
# -   2.) Hugo's Primary cells 2D vs 3D (Clariom D Human Pico) Affymetrix    - #
# ---------------------------------------------------------------------------- #
# Load the processed data
load(file.path(wd, "RData", "HGCC_processedData.RData"))
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
# Save the plot as png
ggsave(file.path(wd, degs_dir, "plots", paste(date, "HGCC-PCA-plot.png", sep = "-")),
       device = "png", HGCC.plots$PCA, 
       bg = "white", dpi = 300, width = 8, height = 6, units = "in")

# Create DEG tables
HGCC.deg <- lapply(HGCC.deg, function(x){
  df = x %>% 
    dplyr::rename(PROBEID = "ID.ID", Symbol = "ID.Symbol", 
                  entrezID = "ID.entrezID", log2FoldChange = "logFC",
                  pvalue = "P.Value", padj = "adj.P.Val")
  
  get_significance(.df = df)
})
# Save tables
sapply(names(HGCC.deg), function(x){
  openxlsx::write.xlsx(HGCC.deg[[x]], 
                       file.path(wd, degs_dir, "tables",
                                 paste(date, "HGCC", x, "DEGs.xlsx", sep = "-")))
})

# Venn diagram from the shared DEGs
HGCC.venn <- list()
HGCC.venn <- get_regions(.list = HGCC.deg[1:3], .names = c("U3017","U3047","U3054"))

(HGCC.plots$VENN <- plot_venn(.data = HGCC.venn$df,
                           .sets = c("U3017","U3047","U3054"),
                           .labels = c("U3054",
                                       "U3047",
                                       "U3017")))
# Save plot
ggsave(file.path(wd, degs_dir, "plots", paste(date, "HGCC-VENN-plot.png", sep = "-")),
       device = "png", HGCC.plots$VENN, 
       bg = "white", dpi = 300, width = 12, height = 8, units = "in")

# Save tables
openxlsx::write.xlsx(HGCC.venn, 
                     file.path(wd, degs_dir, "tables",
                               paste(date, "HGCC-DEGs-venn-regions.xlsx", sep = "-")))

# Create a Vulcano plot for the shared DEGs
(p1 <- plot_vulcan(HGCC.deg$U3017) + ggtitle("U3017"))
(p2 <- plot_vulcan(HGCC.deg$U3047) + ggtitle("U3047"))
(p3 <- plot_vulcan(HGCC.deg$U3054) + ggtitle("U3054"))
(HGCC.plots$vulcano <- cowplot::plot_grid(p1, p2, p3, nrow = 1))
# Save the Vulcano plot
ggsave(file.path(wd, degs_dir, "plots", paste(date, "HGCC-vulcano-plot.png")),
       device = "png", HGCC.plots$vulcano, 
       bg = "white", dpi = 300, width = 20, height = 8, units = "in")

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
sapply(names(HGCC.ORA), function(x){
  openxlsx::write.xlsx(HGCC.ORA[[x]][c(2:3)], 
                       file.path(wd, enrich_dir, "tables",
                                 paste(date, "HGCC", x, "ORA-pathways.xlsx", sep = "-")))
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
                       file.path(wd, enrich_dir, "tables",
                                 paste(date, "HGCC", x, "GSEA-go&hallmark.xlsx", sep = "-")))
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
                       file.path(wd, enrich_dir, "tables",
                                 paste(date, "HGCC", x, "GSEA-pathways.xlsx", sep = "-")))
})

# Save enrichment results
save(HGCC.GO, HGCC.GSEA,
     file = file.path(wd, "RData", "HGCC_enrichment_results.RData"))
# Clear intermediate results
rm(list = c(ls(pattern = "HGCC.*"), "p1","p2","p3","pca_base"))

# Combine enrichment results for CCLD and HGCC cell lines
load(file.path(wd, "RData", "CCLD_enrichment_results.RData"))
load(file.path(wd, "RData", "HGCC_enrichment_results.RData"))

GSEA.object = list(
  U3017 = HGCC.GSEA$U3017$gsea,
  U3047 = HGCC.GSEA$U3047$gsea,
  U3054 = HGCC.GSEA$U3054$gsea,
  CCLD = CCLD.GSEA$gsea
)
GO.object = list(
  U3017 = HGCC.GO$U3017$gsea,
  U3047 = HGCC.GO$U3047$gsea,
  U3054 = HGCC.GO$U3054$gsea,
  CCLD = CCLD.GO$gsea
)
# Visualization of selected enrichment
# 1.) TGF-B: "hsa04350" 
tgfb.ranks = list(
  U3017 = which(HGCC.GSEA$U3017$df$ID == "hsa04350"),
  U3047 = which(HGCC.GSEA$U3047$df$ID == "hsa04350"),
  U3054 = which(HGCC.GSEA$U3054$df$ID == "hsa04350"),
  CCLD = which(CCLD.GSEA$df$ID == "hsa04350"))

col_order = c(4,3,2,1)

HGCC.enrichplot <- list()
HGCC.enrichplot$TGFB$TOTAL_scores  <- do.call(rbind, c(
  lapply(names(tgfb.ranks), function(i){
    gsInfo(tgfb.ranks[[i]], object = GSEA.object[[i]]) %>% 
      dplyr::mutate(source = i)}
  )))

HGCC.enrichplot$TGFB$TOTAL_scores <- HGCC.enrichplot$TGFB$TOTAL_scores %>% 
  dplyr::mutate(
    Description = gsub(x = HGCC.enrichplot$TGFB$TOTAL_scores$Description,
                       pattern = "KEGG_", replacement =  ""),
    source = factor(source, levels = names(GSEA.object)[c(col_order)]))

HGCC.enrichplot$TGFB$TOTAL_table <- getEnrichmentTable(
  .df = rbind(HGCC.GSEA$U3017$df[tgfb.ranks$U3017,] %>% 
                dplyr::mutate(source = "U3017"),
              HGCC.GSEA$U3047$df[tgfb.ranks$U3047,] %>% 
                dplyr::mutate(source = "U3047"),
              HGCC.GSEA$U3054$df[tgfb.ranks$U3054,] %>% 
                dplyr::mutate(source = "U3054"),
              CCLD.GSEA$df[tgfb.ranks$CCLD,] %>% 
                dplyr::mutate(source = "CCLD")),
  .order = c(col_order), .name = "source")

(HGCC.enrichplot$TGFB$TOTAL_plot <- cowplot::plot_grid(
  plotRunningScore(.df = HGCC.enrichplot$TGFB$TOTAL_scores, 
                   .x = "x", .y = "runningScore", 
                   .color = "source", .palette = c("U3017" = "pink",
                                                   "U3047" = "salmon",
                                                   "U3054" = "darkred", 
                                                   "CCLD" = "steelblue"
                   )),
  plotGeneRank(.df = HGCC.enrichplot$TGFB$TOTAL_scores, 
               .x = "x", .facet = "source~.",
               .color = "source", .palette = c("U3017" = "pink",
                                               "U3047" = "salmon",
                                               "U3054" = "darkred",
                                               "CCLD" = "steelblue"
               )),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r",
  margins = c(0.5, 0.5, 0.5, 0.5)))

# Save as svg for publication
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "HGCC-TGFb-GSEA-enrichment-plot.svg", sep = "-")),
       device = "svg", plot = HGCC.enrichplot$TGFB$TOTAL_plot,
       bg = "white", width = 6, height = 4.45, units = "in")
# Save as png for presentation
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "HGCC-TGFb-GSEA-enrichment-plot.png", sep = "-")),
       device = "png", plot = HGCC.enrichplot$TGFB$TOTAL_plot,
       bg = "white", dpi = 300, width = 6, height = 4.45, units = "in")

# 2.) Epithelial-to-mesenchymal transition (EMT): "GO:0001837"
emt.ranks = list(
  U3017 = which(HGCC.GO$U3017$df$ID == "GO:0001837"),
  U3047 = which(HGCC.GO$U3047$df$ID == "GO:0001837"),
  U3054 = which(HGCC.GO$U3054$df$ID == "GO:0001837"),
  CCLD = which(CCLD.GO$df$ID == "GO:0001837"))

HGCC.enrichplot$EMT$TOTAL_scores  <- do.call(rbind, c(
  lapply(names(emt.ranks), function(i){
    gsInfo(emt.ranks[[i]], object = GO.object[[i]]) %>% 
      dplyr::mutate(source = i)}
  )))

HGCC.enrichplot$EMT$TOTAL_scores <- HGCC.enrichplot$EMT$TOTAL_scores %>%
  dplyr::mutate(
    Description = gsub(x = HGCC.enrichplot$EMT$TOTAL_scores$Description,
                       pattern = "GOBP_", replacement =  ""),
    source = factor(source, levels = names(GO.object)[col_order]))

HGCC.enrichplot$EMT$TOTAL_table <- getEnrichmentTable(
  .df = rbind(HGCC.GO$U3017$df[emt.ranks$U3017,] %>% 
                dplyr::mutate(source = "U3017"),
              HGCC.GO$U3047$df[emt.ranks$U3047,] %>% 
                dplyr::mutate(source = "U3047"),
              HGCC.GO$U3054$df[emt.ranks$U3054,] %>% 
                dplyr::mutate(source = "U3054"),
              CCLD.GO$df[emt.ranks$CCLD,] %>% 
                dplyr::mutate(source = "CCLD")),
  .order = col_order, .name = "source")

(HGCC.enrichplot$EMT$TOTAL_plot <- cowplot::plot_grid(
  plotRunningScore(.df = HGCC.enrichplot$EMT$TOTAL_scores, 
                   .x = "x", .y = "runningScore", 
                   .color = "source", .palette = c("U3017" = "pink",
                                                   "U3047" = "salmon",
                                                   "U3054" = "darkred",
                                                   "CCLD" = "steelblue")),
  plotGeneRank(.df = HGCC.enrichplot$EMT$TOTAL_scores, 
               .x = "x", .facet = "source~.",
               .color = "source", .palette = c("U3017" = "pink",
                                               "U3047" = "salmon",
                                               "U3054" = "darkred",
                                               "CCLD" = "steelblue")),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r",
  margins = c(0.5, 0.5, 0.5, 0.5)))

# Save as svg for publication
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "HGCC-EMT-GSEA-enrichment-plot.svg", sep = "-")),
       device = "svg", plot = HGCC.enrichplot$EMT$TOTAL_plot, 
       bg = "white", width = 6, height = 4.45, units = "in")
# Save as png for presentation
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "HGCC-EMT-GSEA-enrichment-plot.png", sep = "-")),
       device = "png", plot = HGCC.enrichplot$EMT$TOTAL_plot, 
       bg = "white", dpi = 300, width = 6, height = 4.45, units = "in")

# 3.) Hypoxia: "M5891"
hypoxia.ranks = list(
  U3017 = which(HGCC.GO$U3017$df$ID == "M5891"),
  U3047 = which(HGCC.GO$U3047$df$ID == "M5891"),
  U3054 = which(HGCC.GO$U3054$df$ID == "M5891"),
  CCLD = which(CCLD.GO$df$ID == "M5891"))

HGCC.enrichplot$HYPOXIA$TOTAL_scores  <- do.call(rbind, c(
  lapply(names(hypoxia.ranks)[1:4], function(i){
    gsInfo(hypoxia.ranks[[i]], object = GO.object[[i]]) %>% 
      dplyr::mutate(source = i)}
  )))

HGCC.enrichplot$HYPOXIA$TOTAL_scores <- HGCC.enrichplot$HYPOXIA$TOTAL_scores %>%
  dplyr::mutate(
    Description = gsub(x = HGCC.enrichplot$HYPOXIA$TOTAL_scores$Description,
                       pattern = "GOBP_", replacement =  ""),
    source = factor(source, levels = names(GO.object)[col_order]))

HGCC.enrichplot$HYPOXIA$TOTAL_table <- getEnrichmentTable(
  .df = rbind(HGCC.GO$U3017$df[hypoxia.ranks$U3017,] %>% 
                dplyr::mutate(source = "U3017"),
              HGCC.GO$U3047$df[hypoxia.ranks$U3047,] %>% 
                dplyr::mutate(source = "U3047"),
              HGCC.GO$U3054$df[hypoxia.ranks$U3054,] %>% 
                dplyr::mutate(source = "U3054"),
              CCLD.GO$df[hypoxia.ranks$CCLD,] %>% 
                dplyr::mutate(source = "CCLD")),
  .order = col_order, .name = "source")

(HGCC.enrichplot$HYPOXIA$TOTAL_plot <- cowplot::plot_grid(
  plotRunningScore(.df = HGCC.enrichplot$HYPOXIA$TOTAL_scores, 
                   .x = "x", .y = "runningScore", 
                   .color = "source", .palette = c("U3017" = "pink",
                                                   "U3047" = "salmon",
                                                   "U3054" = "darkred",
                                                   "CCLD" = "steelblue")),
  plotGeneRank(.df = HGCC.enrichplot$HYPOXIA$TOTAL_scores, 
               .x = "x", .facet = "source~.",
               .color = "source", .palette = c("U3017" = "pink",
                                               "U3047" = "salmon",
                                               "U3054" = "darkred",
                                               "CCLD" = "steelblue")),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r",
  margins = c(0.5, 0.5, 0.5, 0.5)))

# Save as svg for publication
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "HGCC-Hypoxia-GSEA-enrichment-plot.svg", sep = "-")),
       device = "svg", plot = HGCC.enrichplot$HYPOXIA$TOTAL_plot, 
       bg = "white", width = 6, height = 4.45, units = "in")
# Save as png for presentation
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "HGCC-Hypoxia-GSEA-enrichment-plot.png", sep = "-")),
       device = "png", plot = HGCC.enrichplot$HYPOXIA$TOTAL_plot, 
       bg = "white", dpi = 300, width = 6, height = 4.45, units = "in")

# Vulcano visualization of the shared and selected enriched terms
vulcano_pathways <- read.csv("../data/pathways-of-interest.txt", header = T,
                              sep = "\t", stringsAsFactors = T)
vulcano.object <- list(
  U3017 = rbind(HGCC.GSEA$U3017$sig_df[,c("ID","Name", "NES", "p.adjust")],
                HGCC.GO$U3017$sig_df[,c("ID","Name", "NES","p.adjust")]),
  U3047 = rbind(HGCC.GSEA$U3047$sig_df[,c("ID","Name", "NES","p.adjust")],
                HGCC.GO$U3047$sig_df[,c("ID","Name", "NES","p.adjust")]),
  U3054 = rbind(HGCC.GSEA$U3054$sig_df[,c("ID","Name", "NES","p.adjust")],
                HGCC.GO$U3054$sig_df[,c("ID","Name", "NES","p.adjust")]),
  CCLD = rbind(CCLD.GSEA$sig_df[,c("ID","Name", "NES","p.adjust")],
               CCLD.GO$sig_df[,c("ID","Name", "NES","p.adjust")]))

# Cap very low p-values
vulcano.object <- lapply(vulcano.object, function(x){
  x %>% mutate(p.adjust = ifelse(p.adjust < 1e-35, 1e-35, p.adjust))
})

# Count in how many categories each term ID appears
vulcano.ids <- bind_rows(
  lapply(names(vulcano.object), function(name) {
    vulcano.object[[name]] %>%
      select(ID) %>%
      distinct() %>%
      mutate(category = name)
  })
)

# Count appearances per ID
vulcano.ids <- vulcano.ids %>%
  dplyr::group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  dplyr::filter(count >= 3) %>% 
  dplyr::pull(ID)

# Filter vulcano.object to retain only terms present in ≥3 categories
vulcano.object <- lapply(vulcano.object, function(df) {
  df %>% filter(ID %in% vulcano.ids)
})

# 1.) Extracellular matrix organization
HGCC.vulcano <- list()
(HGCC.vulcano$ECM <- lapply(vulcano.object, function(x){
  plotClusters(.df = x, 
               .pathways = vulcano_pathways %>% filter(Category == "ECM"))
}))
# Save individual plots as png
sapply(names(HGCC.vulcano$ECM), function(x){
  ggsave(file.path(wd, enrich_dir, "plots", 
                   paste(date, "HGCC", x, "ECM-GSEA-vulcano-plot.png", sep = "-")),
         device = "png", plot = HGCC.vulcano$ECM[[x]],
         bg = "white", dpi = 300, width = 12, height = 8, units = "in")
})
# Combine the ECM plots into a single figure
(HGCC.vulcano$ECM$total <- cowplot::plot_grid(
  HGCC.vulcano$ECM$U3054 +
    theme(axis.title.x = element_blank()), NULL,
  HGCC.vulcano$ECM$U3047 +
    theme(axis.title = element_blank()), NULL,
  HGCC.vulcano$ECM$U3017 +
    theme(axis.title = element_blank()),
  rel_widths = c(1, 0.15, 1, 0.15, 1),
  byrow = T, nrow = 1, ncol = 5, scale = 1, 
  margins = c(0.5, 0.5, 0.5, 0.5)))

# Save the plot as svg for publication
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "HGCC-ECM-GSEA-vulcano-plot.svg", sep = "-")),
       device = "svg", plot = HGCC.vulcano$ECM$total,
       bg = "white", width = 15, height = 5, units = "in")
# Save the plot as png for presentation
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "HGCC-ECM-GSEA-vulcano-plot.png", sep = "-")),
       device = "png", plot = HGCC.vulcano$ECM$total, 
       bg = "white", dpi = 300, width = 15, height = 5, units = "in")

# 2.) Hypoxia
(HGCC.vulcano$HYPOXIA <- lapply(vulcano.object, function(x){
  plotClusters(.df = x, 
               .pathways = vulcano_pathways %>% filter(Category == "HYPOXIA"))
}))
# Save individual plots
sapply(names(HGCC.vulcano$HYPOXIA), function(x){
  ggsave(file.path(wd, enrich_dir, "plots", 
                   paste(date, "HGCC", x,"Hypoxia-GSEA-vulcano-plot.png",sep = "-")),
         device = "png", plot = HGCC.vulcano$HYPOXIA[[x]], 
         bg = "white", dpi = 300, width = 12, height = 8, units = "in")
})
# Combine the Hypoxia plots into a single figure
(HGCC.vulcano$HYPOXIA$total <- cowplot::plot_grid(
  HGCC.vulcano$HYPOXIA$U3054 +
    theme(axis.title.x = element_blank()), NULL,
  HGCC.vulcano$HYPOXIA$U3047 +
    theme(axis.title = element_blank()), NULL,
  HGCC.vulcano$HYPOXIA$U3017 +
    theme(axis.title = element_blank()),
  rel_widths = c(1, 0.15, 1, 0.15, 1),
  byrow = T, nrow = 1, ncol = 5, scale = 1, 
  margins = c(0.5, 0.5, 0.5, 0.5)))
# Save the plot as svg for publication
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "HGCC-Hypoxia-GSEA-vulcano-plot.svg", sep = "-")),
       device = "svg", plot = HGCC.vulcano$HYPOXIA$total,
       bg = "white", width = 15, height = 5, units = "in")
# Save the plot as png for presentation
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "HGCC-Hypoxia-GSEA-vulcano-plot.png", sep = "-")),
       device = "png", plot = HGCC.vulcano$HYPOXIA$total, 
       bg = "white", dpi = 300, width = 15, height = 5, units = "in")

# 3.) TGF-beta signaling pathway
(HGCC.vulcano$TGFB <- lapply(vulcano.object, function(x){
  plotClusters(.df = x, 
               .pathways = vulcano_pathways %>% filter(Category == "TGFb"))
}))
# Save individual plots
sapply(names(HGCC.vulcano$TGFB), function(x){
  ggsave(file.path(wd, enrich_dir, "plots", 
                   paste(date, "HGCC", x,"TGFb-GSEA-vulcano-plot.png",sep = "-")),
         device = "png", plot = HGCC.vulcano$TGFB[[x]], 
         bg = "white", dpi = 300, width = 12, height = 8, units = "in")
})
# Combine the TGF-b plots into a single figure
(HGCC.vulcano$TGFB$total <- cowplot::plot_grid(
  HGCC.vulcano$TGFB$U3054 +
    theme(axis.title.x = element_blank()), NULL,
  HGCC.vulcano$TGFB$U3047 +
    theme(axis.title = element_blank()), NULL,
  HGCC.vulcano$TGFB$U3017 +
    theme(axis.title = element_blank()),
  rel_widths = c(1, 0.15, 1, 0.15, 1),
  byrow = T, nrow = 1, ncol = 5, scale = 1, 
  margins = c(0.5, 0.5, 0.5, 0.5)))
# Save the plot as svg for publication
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "HGCC-TGFb-GSEA-vulcano-plot.svg", sep = "-")),
       device = "svg", plot = HGCC.vulcano$TGFB$total,
       bg = "white", width = 15, height = 5, units = "in")
# Save the plot as png for presentation
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "HGCC-TGFb-GSEA-vulcano-plot.png", sep = "-")),
       device = "png", plot = HGCC.vulcano$TGFB$total, 
       bg = "white", dpi = 300, width = 15, height = 5, units = "in")


# remove intermediate objects
rm(list = c(ls(pattern = "HGCC.*"), ls(pattern = "CCLD.*"),
            ls(pattern = "*.ranks"), "GO.object", "GSEA.object",
            "vulcano.ids","vulcano.object","col_order"))

# ---------------------------------------------------------------------------- #
# -   3.) U87 Chronic Acidosis AA vs NA & HOX vs NOX (Illumina BeadChip)     - #
# ---------------------------------------------------------------------------- #
# Load the processed data
load(file.path(wd, "RData", "U87_processedData.RData"))
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
pca_base <- prcomp(t(U87.expr[,4:24]), center = TRUE, scale. = TRUE)
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
# Save tables
sapply(names(U87.deg), function(x){
  openxlsx::write.xlsx(U87.deg[[x]], 
                       file.path(wd, degs_dir, "tables",
                                 paste(date, "U87", x, "DEGs.xlsx", sep = "-")))
})

interest_genes_U87 <- c("CSGALNACT1","CSGALNACT2","CHSY1","CHSY3","CHPF","CHPF2",
                        "DCN","BGN","DSE","DSEL","CA9","HIF1A","G0S2","C7orf68")

U87.subset <- list()

U87.subset$logexp <- U87.expr %>%
  dplyr::select(!ENTREZID) %>% 
  dplyr::filter(SYMBOL %in% interest_genes_U87 & 
                  PROBEID %in% U87.deg$`sel_pH647-control_sel`$PROBEID) 

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
                .contrast = c("sel_pH64-control_sel",
                              "acu_pH68-control_acu",
                              "acu_pH64-control_acu",
                              "hypoxia-control_nox"))

tmp <- setNames(tmp, c("sel_pH64","acu_pH68",
                       "acu_pH64","hypoxia"))

U87.subset$df <- lapply(names(tmp), function(x){
  return(tmp[[x]] %>% 
           dplyr::mutate(
             treatment = x,
             SYMBOL = ID.Symbol) %>% 
           dplyr::select(SYMBOL,treatment, logFC))
}) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::full_join(U87.subset$df, ., by = c("SYMBOL", "treatment", "logexp" = "logFC")) %>% 
  dplyr::mutate(logexp = ifelse(is.na(logexp), 0, logexp))

openxlsx::write.xlsx(U87.subset, 
                     file.path(wd, degs_dir, "tables",
                               paste0(date, "-U87-linear-scale-genes-of-interest-expression.xlsx")))

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
sapply(names(U87.ORA), function(x){
  openxlsx::write.xlsx(U87.ORA[[x]][c(2:3)], 
                       file.path(wd, enrich_dir, "tables",
                                 paste(date, "U87", x, "ORA-pathways.xlsx", sep="-")))
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
                       file.path(wd, enrich_dir, "tables",
                                 paste(date, "U87", x,"GSEA-go&hallmark.xlsx", sep = "-")))
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
  openxlsx::write.xlsx(U87.ORA[[x]][c(2:3)], 
                       file.path(wd, degs_dir, "tables",
                                 paste(date, "U87", x, "GSEA-pathways.xlsx", sep = "-")))
})

# Save enrichment results
save(U87.GO, U87.GSEA,
     file = file.path(wd, "RData", "U87_enrichment_results.RData"))
# Clear intermediate results
rm(list = c(ls(pattern = "U87.*"), "interest_genes_U87", "tmp", "pca_base"))

# Combine enrichment results for the U87 cell line
load(file.path(wd, "RData", "U87_enrichment_results.RData"))

GSEA.object = list(
  U87_sel = U87.GSEA$`sel_pH647-control_sel`$gsea,
  U87_acu = U87.GSEA$`acu_pH64-control_acu`$gsea
)
GO.object = list(
  U87_sel = U87.GO$`sel_pH647-control_sel`$gsea,
  U87_acu = U87.GO$`acu_pH64-control_acu`$gsea
)
# Visualization of selected enrichment
# 1.) TGF-B: "hsa04350" 
tgfb.ranks = list(
  "U87_sel" = which(U87.GSEA$`sel_pH647-control_sel`$df$ID == "hsa04350"),
  "U87_acu" = which(U87.GSEA$`acu_pH64-control_acu`$df$ID == "hsa04350")
  )

U87.enrichplot <- list()
U87.enrichplot$TGFB$TOTAL_scores <- do.call(rbind, c(
  lapply(names(tgfb.ranks), function(i){
    gsInfo(tgfb.ranks[[i]], object = GSEA.object[[i]]) %>% 
      dplyr::mutate(source = i)}
  )))
U87.enrichplot$TGFB$TOTAL_scores <- U87.enrichplot$TGFB$TOTAL_scores %>% 
  dplyr::mutate(
    Description = gsub(x = U87.enrichplot$TGFB$TOTAL_scores$Description,
                       pattern = "KEGG_", replacement =  ""),
    source = factor(source, levels = names(GSEA.object)[c(2,1)]))

U87.enrichplot$TGFB$TOTAL_table <- getEnrichmentTable(
  .df = rbind(U87.GSEA$`sel_pH647-control_sel`$df[tgfb.ranks$U87_sel,] %>% 
                dplyr::mutate(source = "U87_sel"),
              U87.GSEA$`acu_pH64-control_acu`$df[tgfb.ranks$U87_acu,] %>% 
                dplyr::mutate(source = "U87_acu")),
  .order = c(2,1), .name = "source")

(U87.enrichplot$TGFB$TOTAL_plot <- cowplot::plot_grid(
  plotRunningScore(.df = U87.enrichplot$TGFB$TOTAL_scores, 
                   .x = "x", .y = "runningScore", 
                   .color = "source", .palette = c("U87_sel" = "#F00000",
                                                   "U87_acu" = "salmon")),
  plotGeneRank(.df = U87.enrichplot$TGFB$TOTAL_scores, 
               .x = "x", .facet = "source~.",
               .color = "source", .palette = c("U87_sel" = "#F00000",
                                               "U87_acu" = "salmon")),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r",
  margins = c(0.5, 0.5, 0.5, 0.5)))
# Save as svg for publication
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "U87-TGFb-GSEA-enrichment-plot.svg", sep = "-")),
       device = "svg", plot = U87.enrichplot$TGFB$TOTAL_plot,
       bg = "white", width = 6, height = 4.45, units = "in")
# Save as png for presentation
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "U87-TGFb-GSEA-enrichment-plot.png", sep = "-")),
       device = "png", plot = U87.enrichplot$TGFB$TOTAL_plot,
       bg = "white", dpi = 300, width = 6, height = 4.45, units = "in")

# 2.) EMT
emt.ranks = list(
  U87_sel = which(U87.GO$`sel_pH647-control_sel`$df$ID == "GO:0001837"),
  U87_acu = which(U87.GO$`acu_pH64-control_acu`$df$ID == "GO:0001837")
  )

U87.enrichplot$EMT$TOTAL_scores <- do.call(rbind, c(
  lapply(names(emt.ranks), function(i){
    gsInfo(emt.ranks[[i]], object = GO.object[[i]]) %>% 
      dplyr::mutate(source = i)}
  )))

U87.enrichplot$EMT$TOTAL_scores <- U87.enrichplot$EMT$TOTAL_scores %>%
  dplyr::mutate(
    Description = gsub(x = U87.enrichplot$EMT$TOTAL_scores$Description,
                       pattern = "GOBP_", replacement =  ""),
    source = factor(source, levels = names(GO.object)[c(2,1)]))

U87.enrichplot$EMT$TOTAL_table <- getEnrichmentTable(
  .df = rbind(U87.GO$`sel_pH647-control_sel`$df[emt.ranks$U87_sel,] %>% 
                dplyr::mutate(source = "U87_sel"),
              U87.GO$`acu_pH64-control_acu`$df[emt.ranks$U87_acu,] %>% 
                dplyr::mutate(source = "U87_acu")),
  .order = c(2,1), .name = "source")

(U87.enrichplot$EMT$TOTAL_plot <- cowplot::plot_grid(
  plotRunningScore(.df = U87.enrichplot$EMT$TOTAL_scores, 
                   .x = "x", .y = "runningScore", 
                   .color = "source", .palette = c("U87_sel" = "#F00000",
                                                   "U87_acu" = "salmon")),
  plotGeneRank(.df = U87.enrichplot$EMT$TOTAL_scores, 
               .x = "x", .facet = "source~.",
               .color = "source", .palette = c("U87_sel" = "#F00000",
                                               "U87_acu" = "salmon")),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r",
  margins = c(0.5, 0.5, 0.5, 0.5)))

# Save as svg for publication
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "U87-EMT-GSEA-enrichment-plot.svg", sep = "-")),
       device = "svg", plot = U87.enrichplot$EMT$TOTAL_plot,
       bg = "white", width = 6, height = 4.45, units = "in")
# Save as png for presentation
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "U87-EMT-GSEA-enrichment-plot.png", sep = "-")),
       device = "png", plot = U87.enrichplot$EMT$TOTAL_plot,
       bg = "white", dpi = 300, width = 6, height = 4.45, units = "in")

# 3.) Hypoxia
hypoxia.ranks = list(
  U87_sel = which(U87.GO$`sel_pH647-control_sel`$df$ID == "M5891"),
  U87_acu = which(U87.GO$`acu_pH64-control_acu`$df$ID == "M5891")
  )

U87.enrichplot$HYPOXIA$TOTAL_scores <- do.call(rbind, c(
  lapply(names(hypoxia.ranks), function(i){
    gsInfo(hypoxia.ranks[[i]], object = GO.object[[i]]) %>% 
      dplyr::mutate(source = i)}
  )))

U87.enrichplot$HYPOXIA$TOTAL_scores <- U87.enrichplot$HYPOXIA$TOTAL_scores %>%
  dplyr::mutate(
    Description = gsub(x = U87.enrichplot$HYPOXIA$TOTAL_scores$Description,
                       pattern = "GOBP_", replacement =  ""),
    source = factor(source, levels = names(GO.object)[c(2,1)]))

U87.enrichplot$HYPOXIA$TOTAL_table <- getEnrichmentTable(
  .df = rbind(U87.GO$`sel_pH647-control_sel`$df[hypoxia.ranks$U87_sel,] %>% 
                dplyr::mutate(source = "U87_CA"),
              U87.GO$`acu_pH64-control_acu`$df[hypoxia.ranks$U87_acu,] %>% 
                dplyr::mutate(source = "U87_AA")),
  .order = c(2,1), .name = "source")

(U87.enrichplot$HYPOXIA$TOTAL_plot <- cowplot::plot_grid(
  plotRunningScore(.df = U87.enrichplot$HYPOXIA$TOTAL_scores, 
                   .x = "x", .y = "runningScore", 
                   .color = "source", .palette = c("U87_sel" = "#F00000",
                                                   "U87_acu" = "salmon")),
  plotGeneRank(.df = U87.enrichplot$HYPOXIA$TOTAL_scores, 
               .x = "x", .facet = "source~.",
               .color = "source", .palette = c("U87_sel" = "#F00000",
                                               "U87_acu" = "salmon")),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r",
  margins = c(0.5, 0.5, 0.5, 0.5)))
# Save as svg for publication
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "U87-Hypoxia-GSEA-enrichment-plot.svg", sep = "-")),
       device = "svg", plot = U87.enrichplot$HYPOXIA$TOTAL_plot, 
       bg = "white", width = 6, height = 4.45, units = "in")
# Save as png for presentation
ggsave(file.path(wd, enrich_dir, "plots", paste(date, "U87-Hypoxia-GSEA-enrichment-plot.png", sep = "-")),
       device = "png", plot = U87.enrichplot$HYPOXIA$TOTAL_plot, 
       bg = "white", dpi = 300, width = 6, height = 4.45, units = "in")

# Vulcano visualization of the shared terms and pathways
U87.vulcano.object <- list(
  U87_CA64 = rbind(U87.GSEA$`sel_pH647-control_sel`$sig_df[,c("ID","Name", "NES","p.adjust")],
                   U87.GO$`sel_pH647-control_sel`$sig_df[,c("ID","Name", "NES","p.adjust")]),
  U87_AA64 = rbind(U87.GSEA$`acu_pH64-control_acu`$sig_df[,c("ID","Name", "NES","p.adjust")],
                   U87.GO$`acu_pH64-control_acu`$sig_df[,c("ID","Name", "NES","p.adjust")]),
  U87_AA68 = rbind(U87.GSEA$`acu_pH68-control_acu`$sig_df[,c("ID","Name", "NES","p.adjust")],
                   U87.GO$`acu_pH68-control_acu`$sig_df[,c("ID","Name", "NES","p.adjust")]),
  U87_Hypoxia = rbind(U87.GSEA$`hypoxia-control_nox`$sig_df[,c("ID","Name", "NES","p.adjust")],
                      U87.GO$`hypoxia-control_nox`$sig_df[,c("ID","Name", "NES","p.adjust")])
)
# Cap very low p-values
U87.vulcano.object <- lapply(U87.vulcano.object, function(x){
  x %>% mutate(p.adjust = ifelse(p.adjust < 1e-35, 1e-35, p.adjust))
})

# Count in how many categories each term ID appears
U87.vulcano.ids <- bind_rows(
  lapply(names(U87.vulcano.object), function(name) {
    U87.vulcano.object[[name]] %>%
      select(ID) %>%
      distinct() %>%
      mutate(category = name)
  })
)

# Count appearances per ID
U87.vulcano.ids <- U87.vulcano.ids %>%
  dplyr::group_by(ID) %>% 
  dplyr::summarise(count = n()) %>%
  dplyr::filter(count >= 3) %>% 
  dplyr::pull(ID)

# Filter vulcano.object to retain only terms present in ≥3 categories
U87.vulcano.object <- lapply(U87.vulcano.object, function(df) {
  df %>% filter(ID %in% U87.vulcano.ids)
})

# 1.) Extracellular matrix organization
U87.vulcano <- list()
(U87.vulcano$ECM <- lapply(U87.vulcano.object, function(x){
  plotClusters(.df = x, 
               .pathways = vulcano_pathways %>% filter(Category == "ECM"))
}))
# Save individual plots
sapply(names(U87.vulcano$ECM), function(x){
  ggsave(file.path(wd, enrich_dir, "plots",
                   paste(date, "U87", x, "ECM-GSEA-Vulcano-plot.png", sep = "-")),
         device = "png", plot = U87.vulcano$ECM[[x]],
         bg = "white", dpi = 300, width = 12, height = 8,units = "in")
})

# 2.) Hypoxia
(U87.vulcano$HYPOXIA <- lapply(U87.vulcano.object, function(x){
  plotClusters(.df = x, 
               .pathways = vulcano_pathways %>% filter(Category == "HYPOXIA"))
}))
# Save individual plots
sapply(names(U87.vulcano$HYPOXIA), function(x){
  ggsave(file.path(wd, enrich_dir, "plots",
            paste(date, "U87", x, "Hypoxia-GSEA-Vulcano-plot.png", sep = "-")),
         device = "png", plot = U87.vulcano$HYPOXIA[[x]],
         bg = "white", dpi = 300, width = 12, height = 8, units = "in")
})

# 3.) TGF-beta signaling pathway
(U87.vulcano$TGFB <- lapply(U87.vulcano.object, function(x){
  plotClusters(.df = x, 
               .pathways = vulcano_pathways %>% filter(Category == "TGFb"))
}))
# Save individual plots
sapply(names(U87.vulcano$TGFB), function(x){
  ggsave(file.path(wd, enrich_dir, "plots",
                   paste(date, "U87", x, "TGFb-GSEA-Vulcano-plot.png", sep = "-")),
         device = "png", plot = U87.vulcano$TGFB[[x]],
         bg = "white", dpi = 300, width = 12, height = 8,units = "in")
})
# remove intermediate objects
rm(list = c(ls(pattern = "U87.*"), ls(pattern = "*.ranks"),
            "GO.object", "GSEA.object"))

# ---------------------------------------------------------------------------- #
# -   4.) PANC1 Chronic Acidosis AA vs NA (Clariom D Human Pico) Affymetrix  - #
# ---------------------------------------------------------------------------- #
load(file.path(wd, "RData", "PANC1_processedData.RData"))

PANC1.deg <- PANC1.deg %>% 
  dplyr::rename(PROBEID = "ID.ID", Symbol = "ID.Symbol",
                entrezID = "ID.entrezID", log2FoldChange = "logFC",
                pvalue = "P.Value", padj = "adj.P.Val") %>% 
  get_significance(.)

# Save tables
openxlsx::write.xlsx(PANC1.deg, file.path(wd, degs_dir, "tables",
                                 paste(date, "PANC1-DEGs.xlsx", sep = "-")))

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

PANC1.subset$df <- tmp[[1]] %>% 
  dplyr::mutate(
    SYMBOL = ID.Symbol) %>% 
  dplyr::select(SYMBOL, logFC) %>% 
  dplyr::bind_rows(.) %>% 
  dplyr::full_join(PANC1.subset$df, ., by = c("SYMBOL")) %>% 
  dplyr::mutate(logFC = ifelse(is.na(logFC), 0, logFC))

openxlsx::write.xlsx(PANC1.subset, 
                     file.path(wd, degs_dir, "tables",
                               paste(date, "PANC1-linear-scale-genes-of-interest-expression.xlsx", sep="-")))

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
                     file.path(wd, enrich_dir, "tables", 
                               paste(date, "PANC1-ORA-pathways.xlsx", sep ="-")))

# GSEA on the GO terms and hallmarks from MSigDB
PANC1.GO <- run_gsea(.geneset = PANC1.genes$background,
                     .terms = terms)
PANC1.GO <- c(PANC1.GO, 
              extract_gsea_results(.gsea = PANC1.GO$gsea, 
                                   .db = terms))
# Save the results
openxlsx::write.xlsx(PANC1.GO[c(2:3)], 
                     file.path(wd, enrich_dir, "tables",
                               paste(date, "PANC1-GSEA-go&hallmark.xlsx", sep ="-")))

# GSEA on the KEGG- and REACTOME pathways from MSigDB
PANC1.GSEA <- run_gsea(.geneset = PANC1.genes$background,
                       .terms = pathways)
PANC1.GSEA <- c(PANC1.GSEA, 
                extract_gsea_results(.gsea = PANC1.GSEA$gsea, 
                                     .db = pathways))
# Save the results
openxlsx::write.xlsx(PANC1.GSEA[c(2:3)],
                     file.path(wd, enrich_dir, "tables",
                               paste(date, "PANC1-GSEA-pathways.xlsx")))

# Save GSEA resutls
save(PANC1.GO, PANC1.GSEA,
     file = file.path(wd, "RData", "PANC1_enrichment_results.RData"))

# Remove intermediate results
rm(list = c(ls(pattern = "PANC.*"), "interest_genes_PANC1", "vulcano_pathways"))

################################################################################
# Clean up the environment                                                     #
################################################################################
# Run garbage collection to free up memory
gc() 
# Clear the environment
rm(list = ls()) 