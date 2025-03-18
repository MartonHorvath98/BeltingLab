################################################################################
# 1.) Calculate Cohen's similarity between all the pathways, terms and hallmark#
################################################################################
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
  dir.create("./data/processed", showWarnings = FALSE)
  write.xlsx(similarity.matrix, file = "./data/processed/term_similarity_matrix.xlsx")
} else {
  load("./RData/term_similarity_matrix.RData")
}



################################################################################
# 2.) Cluster the terms and pathways based on the similarity matrix            #
################################################################################
library(RCy3)
library(igraph)

# -----------        U87 Chronic Acidosis AA vs NA                 ----------- #

# combine GO and pathway enrichment results
U87.combinedGSEA <- list(
  sel_pH64 = rbind(
    dplyr::filter(U87.GSEA$`sel_pH647-control_sel`$sig_df, p.adjust < 0.001),
    dplyr::filter(U87.GO$`sel_pH647-control_sel`$sig_df, p.adjust < 0.001)),
  acu_pH64 = rbind(
    dplyr::filter(U87.GSEA$`acu_pH64-control_acu`$sig_df, p.adjust < 0.001),
    dplyr::filter(U87.GO$`acu_pH64-control_acu`$sig_df, p.adjust < 0.001)))

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
dir.create(file.path("Results","U87",date,"tables"), recursive = T, showWarnings = FALSE)
write.xlsx(U87.cluster$sel_pH64$df, file.path("Results","U87",date,"tables","U87_sel_pH64_GSEA_clusters.xlsx"))
write.xlsx(U87.cluster$acu_pH64$df, file.path("Results","U87",date,"tables","U87_acu_pH64_GSEA_clusters.xlsx"))

# -----------         PANC1 Chronic Acidosis AA vs NA              ----------- #

# Combine GO and pathway enrichment results
PANC1.combinedGSEA <- rbind(
  dplyr::filter(PANC1.GSEA$sig_df, p.adjust < 0.01),
  dplyr::filter(PANC1.GO$sig_df, p.adjust < 0.01))

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
dir.create(file.path("Results","PANC1",date,"tables"), recursive = T, showWarnings = FALSE)
write.xlsx(PANC1.cluster$df, file.path(file.path("Results", "PANC1", date, "tables", "PANC1_GSEA_clusters_p0.01.xlsx")))

################################################################################
# 3. Network visualization of the enriched clusters                            #
################################################################################
# U87 - acute acidosis pH6.4
U87.cluster$acu_pH64$sub_graph <- filter_graph(U87.cluster$acu_pH64$graph, 5)
set.seed(42)
U87.cluster$acu_pH64$layout <- layout_with_fr(U87.cluster$acu_pH64$sub_graph)

(U87.cluster$acu_pH64$plot <- plot_network(.net = U87.cluster$acu_pH64$sub_graph,
                                           .layout = U87.cluster$acu_pH64$layout,
                                           .df = U87.cluster$acu_pH64$df))
# Save plot
ggsave(file.path("Results", "U87", date, "U87_acu_pH64_GSEA_network_NES.png"),
       plot = U87.cluster$acu_pH64$plot, bg = "white",
       width = 20, height = 14, units = "in")
ggsave(file.path("Results", "U87", date, "U87_acu_pH64_GSEA_network_NES.svg"),
       plot = U87.cluster$acu_pH64$plot, bg = "white", device = "svg",
       width = 20, height = 14, units = "in")


# U87 - selective acidosis pH6.4
U87.cluster$sel_pH64$sub_graph <- filter_graph(U87.cluster$sel_pH64$graph, 5)
U87.cluster$sel_pH64$layout <- layout_nicely(U87.cluster$sel_pH64$sub_graph)

(U87.cluster$sel_pH64$plot <- plot_network(.net = U87.cluster$sel_pH64$sub_graph,
                                           .layout = U87.cluster$sel_pH64$layout,
                                           .df = U87.cluster$sel_pH64$df))

# Save plot
ggsave(file.path("Results", "U87", date, "U87_sel_pH64_GSEA_network_NES.png"),
       plot = U87.cluster$sel_pH64$plot, bg = "white",
       width = 20, height = 14, units = "in")
ggsave(file.path("Results", "U87", date, "U87_sel_pH64_GSEA_network_NES.svg"),
       plot = U87.cluster$sel_pH64$plot, bg = "white", device = "svg",
       width = 20, height = 14, units = "in")


# PANC1 - chronic acidosis
PANC1.cluster$sub_graph <- filter_graph(PANC1.cluster$graph, 5)
PANC1.cluster$layout <- layout_nicely(PANC1.cluster$sub_graph)

(PANC1.cluster$plot <- plot_network(.net = PANC1.cluster$sub_graph,
                                    .layout = PANC1.cluster$layout,
                                    .df = PANC1.cluster$df))

# Save plot
ggsave(file.path("Results", "PANC1", date, "PANC1_GSEA_network_NES_p0.01.png"),
       plot = PANC1.cluster$plot, bg = "white",
       width = 20, height = 14, units = "in")
ggsave(file.path("Results", "PANC1", date, "PANC1_GSEA_network_NES_p0.01.svg"),
       plot = PANC1.cluster$plot, bg = "white", device = "svg",
       width = 20, height = 14, units = "in")
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

# ------------------ U87 ------------------ #
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
dir.create(file.path("Results","U87",date,"tables"), recursive = T, showWarnings = FALSE)
write.xlsx(U87.cluster$sel_pH64$interest, file.path("Results","U87",date,"tables","U87_sel_pH64_interest_cluster_genes.xlsx"))
write.xlsx(U87.cluster$acu_pH64$interest, file.path("Results","U87",date,"tables","U87_acu_pH64_interest_cluster_genes.xlsx"))

# U87 - selective acidosis pH6.4
U87.circplot <- list()
U87.circplot$sel_pH64 <- getCircplotData(.cluster = U87.cluster$sel_pH64$df,
                                         .deg = U87.deg$`sel_pH647-control_sel`,
                                         .interest_cluster = interest_cluster, 
                                         .interest_cluster_genes = interest_cluster_genes, 
                                         .palette = cluster_palette)

plotCircplot(.path = file.path("Results", "U87", date, "U87_sel_pH64_GSEA_circosplot.png"),
             .data = U87.circplot$sel_pH64$data.mat,
             .color = U87.circplot$sel_pH64$grid.col,
             .links = U87.circplot$sel_pH64$border.mat,
             .labels = c(interest_cluster, interest_cluster_genes))

# U87 - acute acidosis pH6.4
U87.circplot$acu_pH64 <- getCircplotData(.cluster = U87.cluster$acu_pH64$df,
                                         .deg = U87.deg$`acu_pH64-control_acu`,
                                         .interest_cluster = interest_cluster, 
                                         .interest_cluster_genes = interest_cluster_genes, 
                                         .palette = cluster_palette)

plotCircplot(.path = file.path("Results", "U87", date, "U87_acu_pH64_GSEA_circosplot.png"),
             .data = U87.circplot$acu_pH64$data.mat,
             .color = U87.circplot$acu_pH64$grid.col,
             .links = U87.circplot$acu_pH64$border.mat,
             .labels = c(interest_cluster, interest_cluster_genes))


# ------------------ PANC1 ------------------ #
# Save the contributing genes
PANC1.cluster$interest <- PANC1.cluster$df %>% 
  dplyr::filter(ID %in% names(interest_cluster)) %>% 
  dplyr::select(ID, Name, core_enrichment) %>% 
  dplyr::mutate(core_enrichment = gsub("/",", ", core_enrichment))

# Save results
dir.create(file.path("Results","PANC1",date,"tables"), recursive = T, showWarnings = FALSE)
write.xlsx(PANC1.cluster$interest, file.path("Results","PANC1",date,"tables","PANC1_interest_cluster_genes.xlsx"))

# PANC1 - chronic acidosis
PANC1.circplot <- getCircplotData(.cluster = PANC1.cluster$df,
                                  .deg = PANC1.deg,
                                  .interest_cluster = interest_cluster, 
                                  .interest_cluster_genes = interest_cluster_genes, 
                                  .palette = cluster_palette)

plotCircplot(.path = file.path("Results", "PANC1", date, "PANC1_GSEA_circosplot.png"),
             .data = PANC1.circplot$data.mat,
             .color = PANC1.circplot$grid.col,
             .links = PANC1.circplot$border.mat,
             .labels = c(interest_cluster, interest_cluster_genes))

################################################################################
# 5. Vulcano visualization of the shared terms and pathways                    #
################################################################################

(p1 = plot_vulcan(U87.deg$`sel_pH647-control_sel`, label = F) +
  geom_point(data = subset(U87.deg$`sel_pH647-control_sel`,
                           Symbol %in% interest_cluster_genes),
             aes(x = log2FoldChange, y = -log10(padj)),
             shape = 21, color = "black", fill = "yellow", size = 3, alpha = 0.8) +
  # Add labels for significantly up- or down-regulated genes
  geom_label_repel(data = subset(U87.deg$`sel_pH647-control_sel`,
                          Symbol %in% interest_cluster_genes),
                   aes(x = log2FoldChange, y = -log10(padj), 
                       label = paste(Symbol, " [rank: ", match(subset(U87.deg$`sel_pH647-control_sel`,
                                                                      Symbol %in% interest_cluster_genes)$entrezID,
                                                               names(U87.genes$`sel_pH647-control_sel`$background)), "]", sep = ""),
                                    color = significance)))

(p2 = plot_vulcan(U87.deg$`acu_pH64-control_acu`, label = F) +
  geom_point(data = subset(U87.deg$`acu_pH64-control_acu`,
                           Symbol %in% interest_cluster_genes),
             aes(x = log2FoldChange, y = -log10(padj)),
             shape = 21, color = "black", fill = "yellow", size = 3, alpha = 0.8) +
  # Add labels for significantly up- or down-regulated genes
  geom_label_repel(data = subset(U87.deg$`acu_pH64-control_acu`,
                          Symbol %in% interest_cluster_genes),
                   aes(x = log2FoldChange, y = -log10(padj), 
                       label = paste(Symbol, " [rank: ", match(subset(U87.deg$`acu_pH64-control_acu`,
                                                                      Symbol %in% interest_cluster_genes)$entrezID,
                                                               names(U87.genes$`acu_pH64-control_acu`$background)), "]", sep = ""),
                                    color = significance)))

(p3 = plot_vulcan(PANC1.deg, label = F) +
  geom_point(data = subset(PANC1.deg,
                           Symbol %in% interest_cluster_genes),
             aes(x = log2FoldChange, y = -log10(padj)),
             shape = 21, color = "black", fill = "yellow", size = 3, alpha = 0.8) +
  # Add labels for significantly up- or down-regulated genes
  geom_label_repel(data = subset(PANC1.deg,
                          Symbol %in% interest_cluster_genes),
                   aes(x = log2FoldChange, y = -log10(padj), 
                       label = paste(Symbol, " [rank: ", match(subset(PANC1.deg,
                                                                      Symbol %in% interest_cluster_genes)$entrezID,
                                                               names(PANC1.genes$background)), "]", sep = ""),
                                    color = significance)))

# Save plots
ggsave(file.path("Results", "U87", date, "U87_sel_pH64_GSEA_vulcano.png"),
       plot = p1, bg = "white", width = 14, height = 8, units = "in")
ggsave(file.path("Results", "U87", date, "U87_sel_pH64_GSEA_vulcano.svg"),
       plot = p1, bg = "white", device = "svg", width = 14, height = 8, units = "in")

ggsave(file.path("Results", "U87", date, "U87_acu_pH64_GSEA_vulcano.png"),
       plot = p2, bg = "white", width = 14, height = 8, units = "in")
ggsave(file.path("Results", "U87", date, "U87_acu_pH64_GSEA_vulcano.svg"),
       plot = p2, bg = "white", device = "svg", width = 14, height = 8, units = "in")

ggsave(file.path("Results", "PANC1", date, "PANC1_GSEA_vulcano.png"),
       plot = p3, bg = "white", width = 14, height = 8, units = "in")
ggsave(file.path("Results", "PANC1", date, "PANC1_GSEA_vulcano.svg"),
       plot = p3, bg = "white", device = "svg", width = 14, height = 8, units = "in")
# U87.cluster.vulcano <- list(
#   sel_pH64 = U87.deg$`sel_pH647-control_sel` %>% 
#     dplyr::left_join(linkage, by = c("Symbol" = "node1")) %>%,
#   U3047 = ungroup(TOTAL.clusters) %>%
#     dplyr::filter(grepl("U3047", Regions)) %>% 
#     dplyr::select(ID, Name, U3047.NES, U3047.FDR, Cluster, Type)  %>% 
#     dplyr::rename(NES = U3047.NES, FDR = U3047.FDR),
#   U3054 = ungroup(TOTAL.clusters) %>%
#     dplyr::filter(grepl("U3054", Regions)) %>% 
#     dplyr::select(ID, Name, U3054.NES, U3054.FDR, Cluster, Type)  %>% 
#     dplyr::rename(NES = U3054.NES, FDR = U3054.FDR),
#   CCLD = ungroup(TOTAL.clusters) %>%
#     dplyr::filter(grepl("CCLD", Regions)) %>% 
#     dplyr::select(ID, Name, CCLD.NES, CCLD.FDR, Cluster, Type)  %>% 
#     dplyr::rename(NES = CCLD.NES, FDR = CCLD.FDR))
# 
# interest_pathways <- read.csv("data/pathways-of-interest.txt", header = T,
#                               sep = "\t", stringsAsFactors = T)
# 
# interest_clusters <- dplyr::left_join(interest_pathways, TOTAL.clusters, 
#                                       by = c("ID", "Name")) %>% 
#   dplyr::group_by(Cluster, Category) %>% 
#   dplyr::select(Cluster, Category) %>%
#   dplyr::distinct() 
# 
# cluster_palette = c("89"="#9ECAE1","63"="#4292C6","131"="#0D0887FF",
#                     "146"="#4C02A1FF","186"="#7E03A8FF","2"="#A92395FF",
#                     "90"="#FCBBA1","97"="#CC4678FF","167"="#E56B5DFF",
#                     "159"="#EF3B2C","62"="#F89441FF","164"="#FDC328FF",
#                     "147"="#F0F921FF")
# 
# TOTAL.GSEA.vulcano.plots <- list()
# TOTAL.GSEA.vulcano.plots <- lapply(TOTAL.GSEA.vulcano, function(x){
#   plotClusters(.df = x, .pathways = interest_pathways)
# })
# 
# TOTAL.GSEA.vulcano.plots <- sapply(names(TOTAL.GSEA.vulcano.plots), function(x){
#   ggsave(file.path("Results",paste(x,"GSEA_Vulcano_plot",date,".png",sep = "_")),
#          plot = TOTAL.GSEA.vulcano.plots[[x]], bg = "white",
#          width = 20, height = 14, units = "in")
# })
# 
# table <- TOTAL.GSEA.vulcano$U3017 %>% 
#   dplyr::filter(Name %in% interest_pathways$Name) %>%
#   tibble::column_to_rownames("Name") %>%
#   dplyr::mutate(
#     NES = round(NES, 4),
#     FDR = round(FDR, 4)) %>%
#   dplyr::mutate(FDR = case_when(
#     FDR < 0.001 ~ paste("<0.001", "(***)"),
#     FDR < 0.01 ~ paste(as.character(FDR), "(**)"),
#     FDR < 0.05 ~ paste(as.character(FDR), "(*)"),
#     TRUE ~ as.character(FDR),
#   )) %>% 
#   dplyr::select(!ID)
# 
# table_grob <- tableGrob(
#   CCLD.enrichplot$table, 
#   theme = ttheme_minimal(
#     base_size = 14,
#     core = list(
#       fg_params = list(hjust = 0.5, x = 0.5, col = palette)),
#     rowhead = list(
#       fg_params = list(hjust = 0, x = 0,col = palette[c(5,1,2,3,4)]))
#   )) 
