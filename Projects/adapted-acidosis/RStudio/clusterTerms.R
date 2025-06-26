################################################################################
# Created: 2024 10 25 ; Last Modified: 2025 06 24 ; MH                         #
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

# load MSigDB gene sets
load(file = "./RData/MSigDB_gene_sets.RData")

# Calculate Cohen's similarity between all pathways, terms and hallmark
if (!file.exists("./RData/term_similarity_matrix.RData")){
  # extract a list of gene sets from every reference used:
  total.genesets <- c(split(terms$gene_symbol, # GO terms and MSigDb hallmark sets
                            terms$gs_exact_source),
                      split(pathways$gene_symbol, # KEGG and Reactome pathways
                            pathways$gs_exact_source))
  # calculate the similarity matrix
  similarity.matrix <- cohen_kappa(total.genesets)
  # save the similarity matrix as RData
  save(similarity.matrix, file = "./RData/term_similarity_matrix.RData")
  # save file
  dir.create("../data/processed/similarity_matrix", showWarnings = FALSE)
  write.xlsx(similarity.matrix, ove
             file = "../data/processed/similarity_matrix/term_similarity_matrix.xlsx")
} else {
  load("./RData/term_similarity_matrix.RData")
}
# Get current date for results directory
date <- format(Sys.Date(), "%Y-%m-%d")
# Create results directory
cluster_dir <- "Results/Cluster"
# Create tables and plots directories
dir.create(file.path(wd, cluster_dir, "tables"), recursive = T, showWarnings = FALSE)
dir.create(file.path(wd, cluster_dir, "plots"), recursive = T, showWarnings = FALSE)

################################################################################
#     Cluster the terms and pathways based on the similarity matrix            #
################################################################################
# ---------------------------------------------------------------------------- #
# -                  U87 Chronic Acidosis AA vs NA                           - #
# ---------------------------------------------------------------------------- #
# Load the enrihcment results
load(file = "./RData/U87_enrichment_results.RData")
load(file = "./RData/U87_processedData.RData")
# Combine GO and pathway enrichment results
U87.combinedGSEA <- list(
  sel_pH64 = rbind(
    dplyr::filter(U87.GSEA$`sel_pH647-control_sel`$sig_df, p.adjust < 0.001),
    dplyr::filter(U87.GO$`sel_pH647-control_sel`$sig_df, p.adjust < 0.001)),
  acu_pH64 = rbind(
    dplyr::filter(U87.GSEA$`acu_pH64-control_acu`$sig_df, p.adjust < 0.001),
    dplyr::filter(U87.GO$`acu_pH64-control_acu`$sig_df, p.adjust < 0.001)))
# Preprocess combined GSEA results
U87.combinedGSEA <- lapply(U87.combinedGSEA, function(x){
  x %>%
    dplyr::rowwise(.) %>% 
    # Calculate background ratio
    dplyr::mutate(
      geneRatio = length(unlist(strsplit(core_enrichment,"\\/")))/setSize
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::relocate(c("ID", "Name", "setSize", "geneRatio", "NES", "p.adjust", "core_enrichment"),
                    .before = everything())
})

# Get the clusters
U87.cluster <- lapply(U87.combinedGSEA, get_cluster, similarity.matrix, .threshold = 0.25)

# Get representative terms
U87.cluster$sel_pH64$df <- get_cluster_representative(.cluster = U87.cluster$sel_pH64$df,
                                                      .degs = U87.deg$`sel_pH647-control_sel`)
V(U87.cluster$sel_pH64$graph)$Representative <- U87.cluster$sel_pH64$df$Representative
V(U87.cluster$sel_pH64$graph)$Description <- U87.cluster$sel_pH64$df$Name

U87.cluster$acu_pH64$df <- get_cluster_representative(.cluster = U87.cluster$acu_pH64$df,
                                                      .degs = U87.deg$`acu_pH64-control_acu`)
V(U87.cluster$acu_pH64$graph)$Representative <- U87.cluster$acu_pH64$df$Representative
V(U87.cluster$acu_pH64$graph)$Description <- U87.cluster$acu_pH64$df$Name

# Save results
write.xlsx(U87.cluster$sel_pH64$df, 
           file.path(wd, cluster_dir, "tables",
                     paste(date, "U87-sel_pH64-GSEA-clusters.xlsx", sep = "-")))
write.xlsx(U87.cluster$acu_pH64$df,
           file.path(wd, cluster_dir, "tables",
                     paste(date, "U87_acu_pH64_GSEA_clusters.xlsx", sep = "-")))

# ---------------------------------------------------------------------------- #
# -                  U87 Normoxia vs Hypoxia                                 - #
# ---------------------------------------------------------------------------- #
# combine GO and pathway enrichment results
U87.combinedOX <- rbind(
  dplyr::filter(U87.GSEA$`hypoxia-control_nox`$sig_df, p.adjust < 0.001),
    dplyr::filter(U87.GO$`hypoxia-control_nox`$sig_df, p.adjust < 0.001))
# Preprocess combined GSEA results
U87.combinedOX <- U87.combinedOX %>%
  dplyr::rowwise(.) %>% 
  # Calculate background ratio
  dplyr::mutate(
    geneRatio = length(unlist(strsplit(core_enrichment,"\\/")))/setSize
    ) %>% 
  dplyr::ungroup() %>% 
  dplyr::relocate(c("ID", "Name", "setSize", "geneRatio", "NES", "p.adjust", "core_enrichment"),
                    .before = everything())

# Get the clusters
U87.OX.cluster <- get_cluster(U87.combinedOX, similarity.matrix, .threshold = 0.25)

# Get representative terms
U87.OX.cluster$df <- get_cluster_representative(.cluster = U87.OX.cluster$df,
                                                .degs = U87.deg$`hypoxia-control_nox`)
V(U87.OX.cluster$graph)$Representative <- U87.OX.cluster$df$Representative
V(U87.OX.cluster$graph)$Description <- U87.OX.cluster$df$Name

# Save results
write.xlsx(U87.OX.cluster$df, 
           file.path(wd, cluster_dir, "tables",
                     paste(date,"U87-hypoxia-GSEA-clusters.xlsx", sep = "-")))

# ---------------------------------------------------------------------------- #
# -                  PANC1 Chronic Acidosis AA vs NA                         - #
# ---------------------------------------------------------------------------- #
load(file = "./RData/PANC1_enrichment_results.RData")
load(file = "./RData/PANC1_processedData.RData")
# Combine GO and pathway enrichment results
PANC1.combinedGSEA <- rbind(
  dplyr::filter(PANC1.GSEA$sig_df, p.adjust < 0.01),
  dplyr::filter(PANC1.GO$sig_df, p.adjust < 0.01))
# Preprocess combined GSEA results
PANC1.combinedGSEA <- PANC1.combinedGSEA %>%
  dplyr::rowwise(.) %>% 
  # Calculate background ratio
  dplyr::mutate(
    geneRatio = length(unlist(strsplit(core_enrichment,"\\/")))/setSize
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(ID, Name, setSize, geneRatio, NES, p.adjust, core_enrichment)

# Get the clusters
PANC1.cluster <- get_cluster(PANC1.combinedGSEA, similarity.matrix, .threshold = 0.25)

# Get cluster representatives
PANC1.cluster$df <- get_cluster_representative(.cluster = PANC1.cluster$df,
                                               .degs = PANC1.deg)
V(PANC1.cluster$graph)$Representative <- PANC1.cluster$df$Representative
V(PANC1.cluster$graph)$Description <- PANC1.cluster$df$Name

# Save results
write.xlsx(PANC1.cluster$df, 
           file.path(wd, cluster_dir, "tables",
                     paste(date, "PANC1-GSEA-clusters.xlsx", sep = "-")))

################################################################################
# 3. Network visualization of the enriched clusters                            #
################################################################################
# ---------------------------------------------------------------------------- #
# -                  U87 acute acidosis pH 6.4                               - #
# ---------------------------------------------------------------------------- #
U87.cluster$acu_pH64$sub_graph <- filter_graph(U87.cluster$acu_pH64$graph, 5)
set.seed(42)
U87.cluster$acu_pH64$layout <- layout_with_fr(U87.cluster$acu_pH64$sub_graph)

(U87.cluster$acu_pH64$plot <- plot_network(.net = U87.cluster$acu_pH64$sub_graph,
                                           .layout = U87.cluster$acu_pH64$layout,
                                           .df = U87.cluster$acu_pH64$df))

# Save plot as svg for publication
ggsave(file.path(wd, cluster_dir, "plots",
                 paste(date, "U87-acu_pH64-GSEA-network-NES.svg", sep = "-")),
       device = "svg", plot = U87.cluster$acu_pH64$plot,
       bg = "white", width = 20, height = 14, units = "in")
# Save plot as png for presentation
ggsave(file.path(wd, cluster_dir, "plots",
                 paste(date, "U87-acu_pH64-GSEA-network-NES.png", sep = "-")),
       device = "png", plot = U87.cluster$acu_pH64$plot,
       bg = "white", dpi = 300, width = 20, height = 14, units = "in")

# ---------------------------------------------------------------------------- #
# -                  U87 selective acidosis pH 6.4                           - #
# ---------------------------------------------------------------------------- #
U87.cluster$sel_pH64$sub_graph <- filter_graph(U87.cluster$sel_pH64$graph, 5)
U87.cluster$sel_pH64$layout <- layout_nicely(U87.cluster$sel_pH64$sub_graph)

(U87.cluster$sel_pH64$plot <- plot_network(.net = U87.cluster$sel_pH64$sub_graph,
                                           .layout = U87.cluster$sel_pH64$layout,
                                           .df = U87.cluster$sel_pH64$df))

# Save plot as svg for publication
ggsave(file.path(wd, cluster_dir, "plots",
                 paste(date, "U87-sel_pH64-GSEA-network-NES.svg", sep = "-")),
       device = "svg", plot = U87.cluster$sel_pH64$plot,
       bg = "white", width = 20, height = 14, units = "in")
# Save plot as png for presentation
ggsave(file.path(wd, cluster_dir, "plots",
                 paste(date, "U87-sel_pH64-GSEA-network-NES.png", sep = "-")),
       device = "png", plot = U87.cluster$sel_pH64$plot,
       bg = "white", dpi = 300, width = 20, height = 14, units = "in")

# ---------------------------------------------------------------------------- #
# -                  U87 Hypoxia                                             - #
# ---------------------------------------------------------------------------- #
U87.OX.cluster$sub_graph <- filter_graph(U87.OX.cluster$graph, 5)
set.seed(42)
U87.OX.cluster$layout <- layout_with_fr(U87.OX.cluster$sub_graph)

(U87.OX.cluster$plot <- plot_network(.net = U87.OX.cluster$sub_graph,
                                           .layout = U87.OX.cluster$layout,
                                           .df = U87.OX.cluster$df))
# Save plot as svg for publication
ggsave(file.path(wd, cluster_dir, "plots",
                 paste(date, "U87-hypoxia-GSEA-network-NES.svg", sep = "-")),
       device = "svg", plot = U87.OX.cluster$plot,
       bg = "white", width = 20, height = 14, units = "in")
# Save plot as png for presentation
ggsave(file.path(wd, cluster_dir, "plots",
                 paste(date, "U87-hypoxia-GSEA-network-NES.png", sep = "-")),
       device = "png", plot = U87.OX.cluster$plot,
       bg = "white", dpi = 300, width = 20, height = 14, units = "in")

# ---------------------------------------------------------------------------- #
# -                  PANC1 selective acidosis pH 6.4                         - #
# ---------------------------------------------------------------------------- #
PANC1.cluster$sub_graph <- filter_graph(PANC1.cluster$graph, 5)
PANC1.cluster$layout <- layout_nicely(PANC1.cluster$sub_graph)

(PANC1.cluster$plot <- plot_network(.net = PANC1.cluster$sub_graph,
                                    .layout = PANC1.cluster$layout,
                                    .df = PANC1.cluster$df))

# Save plot as svg for publication
ggsave(file.path(wd, cluster_dir, "plots",
                 paste(date, "PANC1-GSEA-network-NES.svg", sep = "-")),
       device = "svg", plot = PANC1.cluster$plot,
       bg = "white", width = 20, height = 14, units = "in")
# Save plot as png for presentation
ggsave(file.path(wd, cluster_dir, "plots",
                 paste(date, "PANC1-GSEA-network-NES.png", sep = "-")),
       device = "png", plot = PANC1.cluster$plot,
       bg = "white", dpi = 300, width = 20, height = 14, units = "in")

################################################################################
# 4. Circosplot visualization of selected pathways and genes                   #
################################################################################
interest_cluster <- c("R-HSA-1630316" = "GLYCOSAMINOGLYCAN_METABOLISM", 
                      "R-HSA-1793185" = "CHONDROITIN_SULFATE_DERMATAN_SULFATE_METABOLISM", 
                      "GO:0006029" = "PROTEOGLYCAN_METABOLIC_PROCESS",
                      "R-HSA-1474244" = "EXTRACELLULAR_MATRIX_ORGANIZATION",
                      "GO:0045229" = "EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION",
                      "R-HSA-3000178" = "ECM_PROTEOGLYCANS",
                      "M5930" = "EPITHELIAL_MESENCHYMAL_TRANSITION")

cluster_palette <- c("EXTRACELLULAR_MATRIX_ORGANIZATION" = "purple",
                     "PROTEOGLYCAN_METABOLIC_PROCESS" = "green",
                     "GLYCOSAMINOGLYCAN_METABOLISM" = "green",
                     "CHONDROITIN_SULFATE_DERMATAN_SULFATE_METABOLISM" = "green",
                     "EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION" = "purple",
                     "ECM_PROTEOGLYCANS" = "purple",
                     "EPITHELIAL_MESENCHYMAL_TRANSITION" = "purple")

interest_cluster_genes <- c("CA9", "BGN", "DCN", "CSGALNACT1", "SRGN", "CHSY1", "CD44", "CHPF","XYLT1","C7orf68", "G0S2")

# ---------------------------------------------------------------------------- #
# -                  U87 selective- and acute acidosis                       - #
# ---------------------------------------------------------------------------- #
# Save the contributing genes
U87.cluster$sel_pH64$interest <- U87.cluster$sel_pH64$df %>% 
  dplyr::filter(ID %in% names(interest_cluster)) %>% 
  dplyr::select(ID, Name, core_enrichment) %>% 
  dplyr::mutate(core_enrichment = gsub("/",", ", core_enrichment))

U87.cluster$acu_pH64$interest <- U87.cluster$acu_pH64$df %>% 
  dplyr::filter(ID %in% names(interest_cluster)) %>% 
  dplyr::select(ID, Name, core_enrichment) %>% 
  dplyr::mutate(core_enrichment = gsub("/",", ", core_enrichment))

# Save results
write.xlsx(U87.cluster$acu_pH64$interest, 
           file.path(wd, cluster_dir,"tables",
                     paste(date, "U87-acu_pH64-interest_cluster_genes.xlsx", sep = "-")))
write.xlsx(U87.cluster$sel_pH64$interest, 
           file.path(wd, cluster_dir,"tables",
                     paste(date, "U87-sel_pH64-interest_cluster_genes.xlsx", sep = "-")))

# Create circosplot visualization
U87.circplot <- list()
U87.circplot$sel_pH64 <- getCircplotData(.cluster = U87.cluster$sel_pH64$df,
                                         .deg = U87.deg$`sel_pH647-control_sel`,
                                         .interest_cluster = interest_cluster, 
                                         .interest_cluster_genes = interest_cluster_genes, 
                                         .palette = cluster_palette)

plotCircplot(.path = file.path(wd, cluster_dir,"plots",
                               paste(date, "U87-sel_pH64-GSEA-circosplot.png", sep = "-")),
             .data = U87.circplot$sel_pH64$data.mat,
             .color = U87.circplot$sel_pH64$grid.col,
             .links = U87.circplot$sel_pH64$border.mat,
             .labels = c(interest_cluster, interest_cluster_genes))

U87.circplot$acu_pH64 <- getCircplotData(.cluster = U87.cluster$acu_pH64$df,
                                         .deg = U87.deg$`acu_pH64-control_acu`,
                                         .interest_cluster = interest_cluster, 
                                         .interest_cluster_genes = interest_cluster_genes, 
                                         .palette = cluster_palette)

plotCircplot(.path = file.path(wd, cluster_dir,"plots",
                               paste(date, "U87-acu_pH64-GSEA-circosplot.png", sep = "_")),
             .data = U87.circplot$acu_pH64$data.mat,
             .color = U87.circplot$acu_pH64$grid.col,
             .links = U87.circplot$acu_pH64$border.mat,
             .labels = c(interest_cluster, interest_cluster_genes))

# ---------------------------------------------------------------------------- #
# -                  PANC1 selective acidosis pH 6.4                         - #
# ---------------------------------------------------------------------------- #
# Save the contributing genes
PANC1.cluster$interest <- PANC1.cluster$df %>% 
  dplyr::filter(ID %in% names(interest_cluster)) %>% 
  dplyr::select(ID, Name, core_enrichment) %>% 
  dplyr::mutate(core_enrichment = gsub("/",", ", core_enrichment))

# Save results
write.xlsx(PANC1.cluster$interest, 
           file.path(wd, cluster_dir,"tables",
                     paste(date, "PANC1-interest-cluster-genes.xlsx", sep = "-")))

# Create circosplot visualization
PANC1.circplot <- getCircplotData(.cluster = PANC1.cluster$df,
                                  .deg = PANC1.deg,
                                  .interest_cluster = interest_cluster, 
                                  .interest_cluster_genes = interest_cluster_genes, 
                                  .palette = cluster_palette)

plotCircplot(.path = file.path(wd, cluster_dir,"plots",
                               paste(date, "PANC1-GSEA-circosplot.png", sep = "-")),
             .data = PANC1.circplot$data.mat,
             .color = PANC1.circplot$grid.col,
             .links = PANC1.circplot$border.mat,
             .labels = c(interest_cluster, interest_cluster_genes))

################################################################################
# 5. Vulcano visualization of the shared terms and pathways                    #
################################################################################
interest_cluster_genes <- c("CA9", "SRGN", "DCN", "CSGALNACT1", "C7orf68", "BGN", "CHPF", "G0S2")
# Limit the fold change values to a maximum of 6 and minimum of -6 for visualization
U87.cluster.vulcan.df <- U87.deg$`sel_pH647-control_sel` %>% 
  dplyr::rename(
    Symbol = ID.Symbol, 
    log2FoldChange = logFC, 
    padj = adj.P.Val, 
    pvalue = P.Value
  ) %>%
  dplyr::mutate(
    log2FoldChange = ifelse(log2FoldChange > 6, 6,
           ifelse(log2FoldChange < -6, -6, log2FoldChange))
  ) %>% get_significance(.)
  
# Create the vulcano plot of selected genes of interest
(U87.cluster.vulcan = plot_vulcan(U87.cluster.vulcan.df, label = F) +
  geom_point(data = subset(U87.cluster.vulcan.df,
                           Symbol %in% interest_cluster_genes),
             aes(x = log2FoldChange, y = -log10(padj)),
             shape = 21, color = "black", fill = "black", size = 3, alpha = 0.8) +
  scale_x_continuous(limits = c(-6, 6), breaks = seq(-4, 4, 4),
                     expand = expansion(0.01)) +  
  theme(panel.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        panel.grid.major = element_line(color = "grey", linewidth = 0.1),
        panel.grid.minor = element_blank()))

# Save plots
ggsave(file.path(wd, cluster_dir,"plots",
                 paste(date, "U87-sel_pH64-GSEA-vulcano.png")),
       plot = U87.cluster.vulcan, bg = "white", device = "png", dpi = 300,
       width = 8, height = 6, units = "in")

################################################################################
# Clean up the environment                                                     #
################################################################################
# Run garbage collection to free up memory
gc() 
# Clear the environment
rm(list = ls()) 
