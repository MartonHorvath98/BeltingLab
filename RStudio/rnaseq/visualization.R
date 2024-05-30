# ---------------------------------------------------------------------------- #
# -         Differential gene expression - Supplementary figures             - #
# ---------------------------------------------------------------------------- #
###############
# 1. PCA plot #
###############
(pca <- make_pca(deseq.invitro$rld, group = "condition", 
                         labs = c("Tumor tissue",# = "physioxia.Tumour",
                                  "3D spheroid (physioxia)", # = "physioxia.3D",
                                  "Cell culture (2D normoxia)", # = "normoxia.2D",
                                  "Cell culture (2D hypoxia)"), # = "hypoxia.2D"), 
                         cols = c("normoxia.2D" = "steelblue",
                                  "hypoxia.2D" = "salmon",
                                  "physioxia.3D" = "red",
                                  "physioxia.Tumour" = "darkred")))
# Save the plot 
## as .svg file for publication
ggsave(file.path(date, plots_dir,"svg","pca.svg"), plot = pca,
       width = 12, height = 12, units = 'in', device = "svg")
## as .png file for presentation
ggsave(file.path(date, plots_dir,"png","pca.png"), plot = pca,
       width = 12, height = 12, units = 'in', device = "png")
####################
# 2. Vulcano plots #
####################
# ------- In vitro results: using 2D normoxia as baseline expression --------- #
## 2.1 2D hypoxia vs 2D normoxia
(results.2DH_vs_2DN$plot <- make_vulcanplot(
  total = results.2DH_vs_2DN$df, 
  sig = results.2DH_vs_2DN$sig_df))
## 2.2 3D organoids vs 2D normoxia
(results.3D_vs_2DN$plot <- make_vulcanplot(
  total = results.3D_vs_2DN$df, 
  sig = results.3D_vs_2DN$sig_df))
## 2.3 Tumor samples vs 2D normoxia
(results.Tumor_vs_2DN$plot <- make_vulcanplot(
  total = results.Tumor_vs_2DN$df, 
  sig = results.Tumor_vs_2DN$sig_df))

## Merge the results with cowplot
(vulcanplots.invitro <- merge_vulcanplots(
  p1 = results.2DH_vs_2DN$plot, # 2D hypoxia vs 2D normoxia
  p2 = results.3D_vs_2DN$plot, # 3D organoids vs 2D normoxia
  p3 = results.Tumor_vs_2DN$plot, # Tumor samples vs 2D normoxia
  titles = c("2D hypoxia vs 2D normoxia", 
             "3D organoids vs 2D normoxia", 
             "Tumor samples vs 2D normoxia")
))

# Save the plot
## as .svg file for publication
ggsave(file.path(date, plots_dir,"svg","invitro_total_vulcanplot.svg"), 
       plot = vulcanplots.invitro, width = 20, height = 8, units = 'in', device = "svg")
## as .png file for presentation
ggsave(file.path(date, plots_dir,"png","invitro_total_vulcanplot.png"), 
       plot = vulcanplots.invitro, width = 20, height = 8, units = 'in', device = "png")

# -------- Total results: using tumor samples as baseline expression --------- #
## 2.4 2D normoxia vs Tumor
(results.2DN_vs_Tumor$plot <- make_vulcanplot(
  total = results.2DN_vs_Tumor$df, 
  sig = results.2DN_vs_Tumor$sig_df))
## 2.5 2D hypoxia vs Tumor
(results.2DH_vs_Tumor$plot <- make_vulcanplot(
  total = results.2DH_vs_Tumor$df, 
  sig = results.2DH_vs_Tumor$sig_df))
## 2.6 3D organoids vs Tumor
(results.3D_vs_Tumor$plot <- make_vulcanplot(
  total = results.3D_vs_Tumor$df, 
  sig = results.3D_vs_Tumor$sig_df))

## Merge the results with cowplot
(vulcanplots.total <- merge_vulcanplots(
  p1 = results.2DN_vs_Tumor$plot, # 2D normoxia vs Tumor
  p2 = results.2DH_vs_Tumor$plot, # 2D hypoxia vs Tumor
  p3 = results.3D_vs_Tumor$plot, # 3D organoids vs Tumor
  titles = c("2D normoxia vs Tumor", 
             "2D hypoxia vs Tumor", 
             "3D organoids vs Tumor")
))

# Save the plot
## as .svg file for publication
ggsave(file.path(date, plots_dir,"svg","total_vulcanplot.svg"), 
       plot = vulcanplots.total, width = 20, height = 8, units = 'in', device = "svg")
## as .png file for presentation
ggsave(file.path(date, plots_dir,"png","total_vulcanplot.png"), 
       plot = vulcanplots.total, width = 20, height = 8, units = 'in', device = "png")

# ---------------------------------------------------------------------------- #
# -         Figure 1 ) Genome-wide expression module analysis                - #
# ---------------------------------------------------------------------------- #
# 1.) Create heat maps of the normalized expression data
global.heatmap <- make_heatmap(deseq = deseq$rld,
                               expr =  expr.matrix,
                               dend = wgcna$net$dendrograms[[1]], 
                               module = wgcna$modules, 
                               filter = wgcna$net$blockGenes[[1]], 
                               coldata = coldata)

# Save the heatmap
svg(file.path(date, plots_dir,"total_heatmap.svg"),
    width = 18, height = 10)
global.heatmap$plot
dev.off()

# 2.) Module-trait relationships
row_split <- setNames(c("stable","stable", "3D low", "patient", "Tumor low", 
                        "Tumor high", "stable", "stable", "patient", "stable"),
                      row.names(ME_corr))
row_split <- factor(row_split, levels = c("3D low", "Tumor high", "Tumor low", 
                                          "patient","stable"))

col_split <- setNames(c(rep("Condition", 4),
                        rep("Patient", 2)),
                      names(ME_corr))
col_split <- factor(col_split, levels = c("Condition", "Patient"))

(corr_plot <- ComplexHeatmap::Heatmap(
  ME_corr, name = "Pearson's corr.",
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  col = colorRamp2(c(min(ME_corr), 
                     0,
                     max(ME_corr)),
                   c("blue", "white", "red")),
  row_split = row_split,
  column_split = col_split,
  show_row_dend = F, show_row_names = F,
  cluster_row_slices =  T, row_title_rot = 0,
  cluster_columns = F,
  column_names_rot = 45,
  bottom_annotation = HeatmapAnnotation(
    df = data.frame(Condition = colnames(MM)),
    col = list(Condition = c("normoxia.2D" = "steelblue",
                             "hypoxia.2D" = "salmon",
                             "physioxia.3D" = "red",
                             "physioxia.Tumour" = "darkred",
                             "593" = "navy",
                             "673" = "khaki")),
    show_legend = TRUE, show_annotation_name = F
  ),
  left_annotation = rowAnnotation(
    df = data.frame(Module = colnames(MEs)),
    col = list(Module = setNames(colnames(MEs),colnames(MEs))),
    show_legend = TRUE, show_annotation_name = F
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (signif(ME_pval[i, j], 2) < 0.001) {
      grid.text("***", x, y, gp = gpar(fontsize = 12, "bold", col = "black"))
    } else if (signif(ME_pval[i, j], 2) < 0.01) {
      grid.text("**", x, y, gp = gpar(fontsize = 12, "bold", col = "black"))
    } else if (signif(ME_pval[i, j], 2) < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 12, "bold", col = "black"))
    } else {
      grid.text("", x, y, gp = gpar(fontsize = 12, "bold", col = "black"))
    }
    grid.rect(x, y, width, height, gp = gpar(color = "#333333", fill = "transparent",
                                             alpha = 0.5))
  }
))

# Save the results
## as .svg for publication
svg(file.path(date, plots_dir,"svg", "Module_trait_correlation.svg"),
    width = 10, height = 12)
draw(corr_plot)
dev.off()

## Module expression plot
(AE_plot <- ggplot(AE, aes(x = condition, y = averageExpr, fill = cluster)) +
    facet_grid(cluster~., scales = "free_y") +
    geom_violin(position = position_dodge(0.9), # separate by group
                scale = "width", trim = F, # shape configuration
                linewidth = 1, # shape settings
                show.legend = F) + # legend settings
    scale_fill_manual(values = c("Tumor high" = "red",
                                 "Tumor low" = "yellow",
                                 "3D low" = "grey")) +
    # add mean and 1-times the standard deviation statistics to the 
    # violin plots
    stat_summary(aes(group = cluster), # separate by group
                 position = position_dodge(0.9),
                 fun = "mean", geom = "point", # function
                 color = "black", size = 2, # shape settings
                 show.legend = F) + # legend settings
    stat_summary(aes(group = cluster), # separate by group
                 position = position_dodge(0.9), 
                 fun.data="mean_sdl", fun.args = list(mult=1), # function 
                 geom = "errorbar", 
                 color = "black", width = 0.15, size = 1, # shape settings
                 show.legend = F) + # legend settings
    theme_pubclean() + 
    theme(strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(size = 14),
          panel.border = element_rect(color = "black", fill = NA),
          panel.spacing = unit(1, "lines"),
          axis.title = element_blank(),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1)))

# Save the results
## as .svg for publication
ggsave(file.path(date, plots_dir,"svg", "Module_gene_expression.svg"),
       AE_plot, width = 7, height = 12, units = "in")


# ---------------------------------------------------------------------------- #
# -       Figure 2 ) Activation patterns by co-expressed gene modules        - #
# ---------------------------------------------------------------------------- #
# 1.) Compound GO plots
(COMPPLOT.full_modul <- compound_GOplot(COMP.full_modul) + 
    facet_grid(module~category, scales = "free", space = "free_x",
               labeller = as_labeller(c("BP"="Biological processes",
                                        "MF"="Molecular functions",
                                        "CC"="Cell components",
                                        "KEGG"="KEGG pathways",
                                        "3D low"="3D low",
                                        "Tumor high"="Tumor high",
                                        "Tumor low"="Tumor low"))) + 
    theme(panel.spacing = unit(1, "lines")))

# Save the plot
## as .svg file for publication
ggsave(file.path(date, plots_dir,"svg","compound_GOplot_full_modul.svg"),
       COMPPLOT.full_modul, "svg", width = 24, height = 24, units = 'in')
## as .png file for presentation
ggsave(file.path(date, plots_dir,"png","compound_GOplot_full_modul.png"),
       COMPPLOT.full_modul, "png", width = 24, height = 24, units = 'in')

# 2.) Deconvolute stromal cell proportions
## Prepare the data for visualization
cell_proportions <- decon$proportions$qprogwc_sig1 %>% 
  tibble::rownames_to_column("Patient") %>% 
  tidyr::pivot_longer(-Patient, 
                      names_to = "Cells",
                      values_to = "Proportions") %>% 
  dplyr::mutate(Proportions = as.numeric(Proportions/100),
                Patient = factor(Patient, 
                                 levels = c("VI.3429.593.tumor.tissue",
                                            "VI.3429.673.tumor.tissue"),
                                 labels = c("Sample_593","Sample_673")),
                Cells = factor(Cells, 
                               levels = c("astrocytes", "microglia", "macrophages",
                                          "dendritic", "monocytes", "T-cells",
                                          "B-cells", "GBM.cells"),
                               labels = c("Astrocytes", "Microglial cells", 
                                          "Macrophages", "Dendritic cells", 
                                          "Monocytes", "T-cells", "B-cells",
                                          "GBM cells")))

# Define the colors for the cell types
cell_colors <- c(
  "Astrocytes" = "#ADD8E6",
  "Microglial cells" = "#00008B",
  "Macrophages" = "#90EE90",
  "Dendritic cells" = "#32CD32",
  "Monocytes" = "#006400",
  "T-cells" = "#FFA07A",
  "B-cells" = "#8B0000"
)
# Create the plot
p <- ggplot(data = cell_proportions,
            aes(y = Patient, x = Proportions, fill = Cells)) +
  geom_bar(position = "fill", stat = "identity", color='black',width=0.9) +
  geom_text(data = cell_proportions %>% dplyr::filter(Proportions != 0),
            aes(label = scales::percent(Proportions)),
            position = position_stack(vjust = 0.5)) +
  scale_x_continuous(labels = scales::percent, expand = c(0,0),
                     position = "top", name = "Estimated cell type proportions (%)") +
  scale_fill_manual(name = "Cell types",
                    values = cell_colors) +
  theme_classic() + 
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

# Create separate legends
glial_cells_legend <- ggplot(
  data = cell_proportions %>% dplyr::filter(Cells %in% c("Astrocytes", "Microglial cells")),
  aes(y = Patient, x = Proportions, fill = Cells)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Glial Cells", values = cell_colors[c("Astrocytes", "Microglial cells")]) +
  theme_void() +
  theme(legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))

monocyte_derived_legend <-  ggplot(
  data = cell_proportions %>% dplyr::filter(Cells %in% c("Macrophages", "Dendritic cells", "Monocytes")),
  aes(y = Patient, x = Proportions, fill = Cells)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Monocyte-derived Cells", values = cell_colors[c("Macrophages", "Dendritic cells", "Monocytes")]) +
  theme_void() +
  theme(legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))

lymphocytes_legend <- ggplot(
  data = cell_proportions %>% dplyr::filter(Cells %in% c("T-cells", "B-cells")),
  aes(y = Patient, x = Proportions, fill = Cells)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Lymphocytes", values = cell_colors[c("T-cells", "B-cells")]) +
  theme_void()  +
  theme(legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12))

# Extract the legends
glial_cells_legend_grob <- get_legend(glial_cells_legend)
monocyte_derived_legend_grob <- get_legend(monocyte_derived_legend)
lymphocytes_legend_grob <- get_legend(lymphocytes_legend)

# Combine the legends
combined_legends <- cowplot::plot_grid(
  glial_cells_legend_grob,
  monocyte_derived_legend_grob,
  lymphocytes_legend_grob,
  ncol = 1,
  align = "hv"
)

## Load PUREE results
puree <- read.table(file.path(date, results_dir, "PUREE_results.txt"),
                    header = T, sep = "\t")
## Create PUREE plot
q <- ggplot(data = puree,
            aes(y = Patient, x = Purity)) +
  geom_bar(stat = "identity", fill = "orange", color='black',width=0.9) +
  geom_text(aes(label = scales::percent(Purity)),
            position = position_stack(vjust = 0.5)) +
  scale_x_continuous(labels = scales::percent,
                     limits = c(0,1), expand = c(0,0),
                     position = "top", name = "Estimated tumor purity (%)") +
  theme_classic() + 
  theme(legend.position = "none",
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

# Combine the base plot with the legends
final_plot <- cowplot::plot_grid( q, p, combined_legends,
                                  ncol = 3, rel_widths = c(2, 2, 1))

# Display the final plot
ggsave(file.path(date, plots_dir, "cell_proportions.svg"),
       final_plot, width = 16, height = 8, units = "in", device = "svg")

# ---------------------------------------------------------------------------- #
# -         Pathway analysis visualization - Supplementary figures           - #
# ---------------------------------------------------------------------------- #
# 1.) GO Semantic similarity visual analysis 
(GOSS.2DN_vs_Tumor <- make_GO_simplot(REDUCED.2DN_vs_Tumor$BP$reduced,
                                      REDUCED.2DN_vs_Tumor$BP$subset))
ggsave(file.path(date, plots_dir,"2DN_vs_Tumor_GO_BP.png"),
       plot = GOSS.2DN_vs_Tumor, width = 16, height = 12, units = 'in')

## 1.2 2D normoxia vs Tumor - Cellular Component
(GOSS.2DN_vs_Tumor <- make_GO_simplot(REDUCED.2DN_vs_Tumor$CC$reduced,
                                      REDUCED.2DN_vs_Tumor$CC$subset))
ggsave(file.path(date, plots_dir,"2DN_vs_Tumor_GO_CC.png"),
       plot = GOSS.2DN_vs_Tumor, width = 16, height = 12, units = 'in')

## 1.3 2D normoxia vs Tumor - Molecular Function
(GOSS.2DN_vs_Tumor <- make_GO_simplot(REDUCED.2DN_vs_Tumor$MF$reduced,
                                      REDUCED.2DN_vs_Tumor$MF$subset))
ggsave(file.path(date, plots_dir,"2DN_vs_Tumor_GO_MF.png"),
       plot = GOSS.2DN_vs_Tumor, width = 16, height = 12, units = 'in')

## 2.1 2D hypoxia vs Tumor - Biological Process
(GOSS.2DH_vs_Tumor <- make_GO_simplot(REDUCED.2DH_vs_Tumor$BP$reduced,
                                      REDUCED.2DH_vs_Tumor$BP$subset))
ggsave(file.path(date, plots_dir,"2DH_vs_Tumor_GO_BP.png"),
       plot = GOSS.2DH_vs_Tumor, width = 16, height = 12, units = 'in')

## 2.2 2D hypoxia vs Tumor - Cellular Component
(GOSS.2DH_vs_Tumor <- make_GO_simplot(REDUCED.2DH_vs_Tumor$CC$reduced,
                                      REDUCED.2DH_vs_Tumor$CC$subset))
ggsave(file.path(date, plots_dir,"2DH_vs_Tumor_GO_CC.png"),
       plot = GOSS.2DH_vs_Tumor, width = 16, height = 12, units = 'in')

## 2.3 2D hypoxia vs Tumor - Molecular Function
(GOSS.2DH_vs_Tumor <- make_GO_simplot(REDUCED.2DH_vs_Tumor$MF$reduced,
                                      REDUCED.2DH_vs_Tumor$MF$subset))
ggsave(file.path(date, plots_dir,"2DH_vs_Tumor_GO_MF.png"),
       plot = GOSS.2DH_vs_Tumor, width = 16, height = 12, units = 'in')

## 3.1 3D organoids vs Tumor - Biological Process
(GOSS.3D_vs_Tumor <- make_GO_simplot(REDUCED.3D_vs_Tumor$BP$reduced,
                                     REDUCED.3D_vs_Tumor$BP$subset))
ggsave(file.path(date, plots_dir,"3D_vs_Tumor_GO_BP.png"),
       plot = GOSS.3D_vs_Tumor, width = 16, height = 12, units = 'in')

## 3.2 3D organoids vs Tumor - Cellular Component
(GOSS.3D_vs_Tumor <- make_GO_simplot(REDUCED.3D_vs_Tumor$CC$reduced,
                                     REDUCED.3D_vs_Tumor$CC$subset))
ggsave(file.path(date, plots_dir,"3D_vs_Tumor_GO_CC.png"),
       plot = GOSS.3D_vs_Tumor, width = 16, height = 12, units = 'in')

## 3.3 3D organoids vs Tumor - Molecular Function
(GOSS.3D_vs_Tumor <- make_GO_simplot(REDUCED.3D_vs_Tumor$MF$reduced,
                                     REDUCED.3D_vs_Tumor$MF$subset))
ggsave(file.path(date, plots_dir,"3D_vs_Tumor_GO_MF.png"),
       plot = GOSS.3D_vs_Tumor, width = 16, height = 12, units = 'in')

###################
# 3. GO SS plots  #
###################
## 4.1 2D hypoxia vs 2D normoxia - Biological Process
(GOSS.2DH_vs_2DN <- make_GO_simplot(REDUCED.2DH_vs_2DN$BP$reduced,
                                   REDUCED.2DH_vs_2DN$BP$subset))
ggsave(file.path(date, plots_dir,"2DH_vs_2DN_GO_BP.png"),
       plot = GOSS.2DH_vs_2DN, width = 16, height = 12, units = 'in')

## 4.2 2D hypoxia vs 2D normoxia - Cellular Component
(GOSS.2DH_vs_2DN <- make_GO_simplot(REDUCED.2DH_vs_2DN$CC$reduced,
                                   REDUCED.2DH_vs_2DN$CC$subset))
ggsave(file.path(date, plots_dir,"2DH_vs_2DN_GO_CC.png"),
       plot = GOSS.2DH_vs_2DN, width = 16, height = 12, units = 'in')

## 4.3 2D hypoxia vs 2D normoxia - Molecular Function
(GOSS.2DH_vs_2DN <- make_GO_simplot(REDUCED.2DH_vs_2DN$MF$reduced,
                                   REDUCED.2DH_vs_2DN$MF$subset))
ggsave(file.path(date, plots_dir,"2DH_vs_2DN_GO_MF.png"),
       plot = GOSS.2DH_vs_2DN, width = 16, height = 12, units = 'in')

## 5.1 3D organoids vs 2D normoxia - Biological Process
(GOSS.3D_vs_2DN <- make_GO_simplot(REDUCED.3D_vs_2DN$BP$reduced,
                                   REDUCED.3D_vs_2DN$BP$subset))
ggsave(file.path(date, plots_dir,"3D_vs_2DN_GO_BP.png"),
       plot = GOSS.3D_vs_2DN, width = 16, height = 12, units = 'in')

## 5.2 3D organoids vs 2D normoxia - Cellular Component
(GOSS.3D_vs_2DN <- make_GO_simplot(REDUCED.3D_vs_2DN$CC$reduced,
                                   REDUCED.3D_vs_2DN$CC$subset))
ggsave(file.path(date, plots_dir,"3D_vs_2DN_GO_CC.png"),
       plot = GOSS.3D_vs_2DN, width = 16, height = 12, units = 'in')

## 5.3 3D organoids vs 2D normoxia - Molecular Function
(GOSS.3D_vs_2DN <- make_GO_simplot(REDUCED.3D_vs_2DN$MF$reduced,
                                   REDUCED.3D_vs_2DN$MF$subset))

ggsave(file.path(date, plots_dir,"3D_vs_2DN_GO_MF.png"),
       plot = GOSS.3D_vs_2DN, width = 16, height = 12, units = 'in')

#########################
# 5.) Compound GO plots #
#########################
## 1. 2D normoxia vs Tumor
(COMPPLOT.2DN_vs_Tumor <- (compound_GOplot(COMP.2DN_vs_Tumor)))
ggsave(file.path(date, plots_dir,"2DN_vs_Tumor_compound_GO.png"),
       COMPPLOT.2DN_vs_Tumor, "png", width = 20, height = 12, units = 'in')

## 2. 2D hypoxia vs Tumor
(COMPPLOT.2DH_vs_Tumor <- (compound_GOplot(COMP.2DH_vs_Tumor)))
ggsave(file.path(date, plots_dir,"2DH_vs_Tumor_compound_GO.png"),
       COMPPLOT.2DH_vs_Tumor, "png", width = 20, height = 12, units = 'in')

## 3. 3D organoids vs Tumor
(COMPPLOT.3D_vs_Tumor <- (compound_GOplot(COMP.3D_vs_Tumor)))
ggsave(file.path(date, plots_dir,"3D_vs_Tumor_compound_GO.png"),
       COMPPLOT.3D_vs_Tumor, "png", width = 20, height = 12, units = 'in')

## 4. 2D hypoxia vs 2D normoxia
(COMPPLOT.2DH_vs_2DN <- (compound_GOplot(COMP.2DH_vs_2DN)))
ggsave(file.path(date, plots_dir,"2DH_vs_2DN_compound_GO.png"),
       COMPPLOT.2DH_vs_2DN, "png", width = 20, height = 12, units = 'in')

## 5. 3D organoids vs 2D normoxia
(COMPPLOT.3D_vs_2DN <- (compound_GOplot(COMP.3D_vs_2DN)))
ggsave(file.path(date, plots_dir,"3D_vs_2DN_compound_GO.png"),
       COMPPLOT.3D_vs_2DN, "png", width = 20, height = 12, units = 'in')


# ------------------------------------------------------------------------------
# -         Figure 3.) COMPARING PATIENT 593 AND PATIENT 673                   -
# ------------------------------------------------------------------------------
###################
# 1. Vulcan plots #
###################
# ------- Patient 593: using 2D normoxia as baseline expression --------- #
## 1.1 2D hypoxia vs 2D normoxia
(results.593.2DH_vs_2DN$plot <- make_vulcanplot(
  total = results.593.2DH_vs_2DN$df, 
  sig = results.593.2DH_vs_2DN$sig_df))
## 1.2 3D organoids vs 2D normoxia
(results.593.3D_vs_2DN$plot <- make_vulcanplot(
  total = results.593.3D_vs_2DN$df, 
  sig = results.593.3D_vs_2DN$sig_df))
## 1.3 Tumor samples vs 2D normoxia
(results.593.Tumor_vs_2DN$plot <- make_vulcanplot(
  total = results.593.Tumor_vs_2DN$df, 
  sig = results.593.Tumor_vs_2DN$sig_df))

## Merge the results with cowplot
(vulcanplots.593 <- merge_vulcanplots(
  p1 = results.593.2DH_vs_2DN$plot, # 2D hypoxia vs 2D normoxia
  p2 = results.593.3D_vs_2DN$plot, # 3D organoids vs 2D normoxia
  p3 = results.593.Tumor_vs_2DN$plot, # Tumor samples vs 2D normoxia
  titles = c("2D hypoxia vs 2D normoxia", 
             "3D organoids vs 2D normoxia", 
             "Tumor samples vs 2D normoxia")
))

# Save the plot
## as .svg file for publication
ggsave(file.path(date, plots_dir,"svg","patient_593_vulcanplot.svg"), 
       plot = vulcanplots.593, width = 20, height = 8, units = 'in', device = "svg")
## as .png file for presentation
ggsave(file.path(date, plots_dir,"png","patient_593_vulcanplot.png"), 
       plot = vulcanplots.593, width = 20, height = 8, units = 'in', device = "png")

# ------- Patient 673: using 2D normoxia as baseline expression --------- #
## 2.4 2D hypoxia vs 2D normoxia
(results.673.2DH_vs_2DN$plot <- make_vulcanplot(
  total = results.673.2DH_vs_2DN$df, 
  sig = results.673.2DH_vs_2DN$sig_df))
## 2.5 3D organoids vs 2D normoxia
(results.673.3D_vs_2DN$plot <- make_vulcanplot(
  total = results.673.3D_vs_2DN$df, 
  sig = results.673.3D_vs_2DN$sig_df))
## 2.6 Tumor samples vs 2D normoxia
(results.673.Tumor_vs_2DN$plot <- make_vulcanplot(
  total = results.673.Tumor_vs_2DN$df, 
  sig = results.673.Tumor_vs_2DN$sig_df))

## Merge the results with cowplot
(vulcanplots.673 <- merge_vulcanplots(
  p1 = results.673.2DH_vs_2DN$plot, # 2D hypoxia vs 2D normoxia
  p2 = results.673.3D_vs_2DN$plot, # 3D organoids vs 2D normoxia
  p3 = results.673.Tumor_vs_2DN$plot, # Tumor samples vs 2D normoxia
  titles = c("2D hypoxia vs 2D normoxia", 
             "3D organoids vs 2D normoxia", 
             "Tumor samples vs 2D normoxia")
))

# Save the plot
## as .svg file for publication
ggsave(file.path(date, plots_dir,"svg","patient_673_vulcanplot.svg"), 
       plot = vulcanplots.673, width = 20, height = 8, units = 'in', device = "svg")
## as .png file for presentation
ggsave(file.path(date, plots_dir,"png","patient_673_vulcanplot.png"), 
       plot = vulcanplots.673, width = 20, height = 8, units = 'in', device = "png")

#################
# Venn diagram  #
#################
## Patient 593 
LFC_range <- range(min(c(upset.593$log2FoldChange, upset.673$log2FoldChange)),
                   max(c(upset.593$log2FoldChange, upset.673$log2FoldChange)))

(venn_plot.593 <- ( # make the plot
  ggplot(upset.593) 
  + coord_fixed()
  + theme_void() # remove the background
  + geom_point(aes(x = x, y = y, colour = log2FoldChange), size = 2) # add the genes
  + scale_color_gradient2(
    limits = LFC_range,
    low = "blue", mid = "white", high = "red", midpoint = 0,
    na.value = "grey50", guide = "colourbar", aesthetics = "colour"
  ) # add a color gradient to the points
  + geom_venn_region(upset.593, sets = c("2DN","2DH", "3D"), 
                     alpha = 0.3, show.legend = F) # add the regions
  + geom_venn_circle(upset.593, sets = c("2DN","2DH", "3D"),
                     size = .5) # add the circles
  + scale_fill_venn_mix(upset.593, sets = c("2DN","2DH", "3D"),
                        guide='none', highlight=c("2DN-2DH-3D"),
                        inactive_color='NA') # highlight the 3-way intersection
  + geom_venn_label_set(upset.593, outwards_adjust = 2,
                        sets = c("3D","2DH", "2DN"),
                        fill = alpha("black", .35), 
                        aes(label = c("2DN","2DH", "3D")))) 
  # add the set labels
  + geom_venn_label_region(
    upset.593, sets = c("2DN","2DH", "3D"),
    aes(label=size), outwards_adjust=1.25, position=position_nudge(y=0.2)) 
  # add the number of genes in each region
  + theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ))

# save the plot
ggsave(file.path(date, plots_dir,"png","venn_plot.593.png"), venn_plot.593,
       "png", width = 12, height = 6, units = 'in')

## --------------- Patient 673 -------------------------------------------------
(venn_plot.673 <- ( # make the plot
  ggplot(upset.673) 
  + coord_fixed()
  + theme_void() # remove the background
  + geom_point(aes(x = x, y = y, colour = log2FoldChange), size = 2) # add the genes
  + scale_color_gradient2(
    limits = LFC_range,
    low = "blue", mid = "white", high = "red", midpoint = 0,
    na.value = "grey50", guide = "colourbar", aesthetics = "colour"
  ) # add a color gradient to the points
  + geom_venn_region(upset.673, sets = c("2DN","2DH", "3D"), 
                     alpha = 0.3, show.legend = F) # add the regions
  + geom_venn_circle(upset.673, sets = c("2DN","2DH", "3D"),
                     size = .5) # add the circles
  + scale_fill_venn_mix(upset.673, sets = c("2DN","2DH", "3D"),
                        guide='none', highlight = c("2DN-2DH-3D"),
                        inactive_color='NA') # highlight the 3-way intersection
  + geom_venn_label_set(upset.673, outwards_adjust = 2,
                        sets = c("3D","2DH", "2DN"),
                        fill = alpha("black", .35), 
                        aes(label = c("2DN","2DH", "3D")))) 
 # add the set labels
 + geom_venn_label_region(
   upset.673, sets = c("2DN","2DH", "3D"),
   aes(label=size), outwards_adjust=1.25, position=position_nudge(y=0.2)) 
 # add the number of genes in each region
 + theme(
   legend.title = element_text(size = 12),
   legend.text = element_text(size = 12)
 ))

# Save the plot
ggsave(file.path(date, plots_dir, "venn_plot.673.png"), venn_plot.673,
       "png", width = 12, height = 6, units = 'in')

venn.legend <- get_legend(venn_plot.593)

venn.grid <- cowplot::plot_grid(venn_plot.593 + 
                     theme(legend.position = "none"),
                   venn_plot.673 + 
                     theme(legend.position = "none"),
                   venn.legend, 
                   labels = c("Patient 593","Patient 673"),
                   ncol = 3, rel_widths = c(1,1,0.3))

# Save the plot
ggsave(file.path(date, plots_dir, "png", "venn_plot.png"), venn.grid,
       "png", width = 20, height = 8, units = 'in')
ggsave(file.path(date, plots_dir, "svg", "venn_plot.svg"), venn.grid,
       "svg", width = 20, height = 8, units = 'in')


################################################
# ORA visualization                            #
################################################
# Create the ORA plot
ORA.593.compare <- rbind(ORA.593.2DN_vs_Tumor %>% dplyr::mutate(cluster = "2DN"),
                         ORA.593.2DH_vs_Tumor %>% dplyr::mutate(cluster = "2DH"),
                         ORA.593.3D_vs_Tumor %>% dplyr::mutate(cluster = "3D")) %>%
  dplyr::mutate(cluster = factor(cluster, levels = c("2DN", "2DH", "3D")))

ORA.593.compare.df <- ORA.593.compare %>%
  dplyr::select(-geneID) %>% 
  dplyr::mutate(Description = forcats::fct_reorder(Description, zscore)) %>% 
  group_by(cluster, category) %>%
  dplyr::slice_head(n = 5)

ORA.673.compare <- rbind(ORA.673.2DN_vs_Tumor %>% dplyr::mutate(cluster = "2DN"),
                         ORA.673.2DH_vs_Tumor %>% dplyr::mutate(cluster = "2DH"),
                         ORA.673.3D_vs_Tumor %>% dplyr::mutate(cluster = "3D")) %>%
  dplyr::mutate(cluster = factor(cluster, levels = c("2DN", "2DH", "3D")))

ORA.673.compare.df <- ORA.673.compare %>%
  dplyr::select(-geneID) %>% 
  dplyr::mutate(Description = forcats::fct_reorder(Description, zscore)) %>% 
  group_by(cluster, category) %>%
  dplyr::slice_head(n = 5)


size_limits <- range(min(c(ORA.593.compare$Count, ORA.673.compare$Count)),
                     max(c(ORA.593.compare$Count, ORA.673.compare$Count)))
fill_limits <- range(min(c(ORA.593.compare$zscore, ORA.673.compare$zscore)),
                     max(c(ORA.593.compare$zscore, ORA.673.compare$zscore)))

(ORA.593.plot <- ggplot(ORA.593.compare.df, aes(x = cluster, y = Description, 
                               fill = zscore, size = Count)) +
  geom_point(shape = 21, color = "black") + 
  scale_fill_gradient2(
    limits = fill_limits,
    low = "blue", mid = "white", high = "red", midpoint = 0,
    na.value = "grey50", guide = "colourbar", aesthetics = "fill"
  ) +
  scale_size(limits = size_limits,range = c(2,7)) + 
  guides(size = guide_legend(title = "Gene count", order = 1),
         fill = guide_colorbar(title = "Activation (z-score)", order = 2)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) + 
  facet_grid(category~., drop = T, scales = "free_y") + 
  theme_bw() + 
  theme(legend.position = "right",
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")))

# Save the plot
ggsave(file.path(date, plots_dir,"png", "ORA_593_compare.png"), ORA.593.plot,
       "png", width = 8, height = 14, units = 'in')

(ORA.673.plot <- ggplot(ORA.673.compare.df, aes(x = cluster, y = Description, 
                                               fill = zscore, size = Count)) +
  geom_point(shape = 21, color = "black") + 
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    na.value = "grey50", guide = "colourbar", aesthetics = "fill"
  ) +
  scale_size(range = c(2,7)) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) + 
  facet_grid(category~., drop = T, scales = "free_y") + 
  guides(size = guide_legend(title = "Gene count", order = 1),
         fill = guide_colorbar(title = "Activation (z-score)", order = 2)) +
  theme_bw() + 
  theme(legend.position = "right",
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")))

# Save the plot
ggsave(file.path(date, plots_dir,"png", "ORA_673_compare.png"), ORA.673.plot,
       "png", width = 8, height = 14, units = 'in')

ora.legend <- get_legend(ORA.593.plot)

(ora.grid <- cowplot::plot_grid(ORA.593.plot + 
                                  theme(legend.position = "none"),
                               ORA.673.plot + 
                                  theme(legend.position = "none"),
                                ora.legend, 
                                #labels = c("Patient 593","Patient 673"),
                                ncol = 3, rel_widths = c(1,1,0.3)))

# Save the plot
ggsave(file.path(date, plots_dir,"png", "ORA_compare.png"), ora.grid,
       "png", width = 14, height = 12, units = 'in')
# Save the plot
ggsave(file.path(date, plots_dir,"svg", "ORA_compare.svg"), ora.grid,
       "svg", width = 14, height = 12, units = 'in')


################################################
# 8.) SURFME genes visualization               #
################################################
# ---------------- Total results -----------------------------------------------
# 1. 2D normoxia vs Tumor
surfme.2DN_vs_Tumor <- get_CategoryExpressionPlot(results.2DN_vs_Tumor$df, 
                                                  genes = surfme_diff$`2DN`$gene)
ggsave(file.path(date, plots_dir, "2DN_vs_Tumor_surfme.png"),
       surfme.2DN_vs_Tumor, "png", width = 10, height = 8, units = 'in')

# 2. 2D hypoxia vs Tumor
surfme.2DH_vs_Tumor <- get_CategoryExpressionPlot(results.2DH_vs_Tumor$df, 
                                                  genes = surfme_diff$`2DH`)
ggsave(file.path(date, plots_dir, "2DH_vs_Tumor_surfme.png"),
       surfme.2DH_vs_Tumor, "png", width = 10, height = 8, units = 'in')

# 3. 3D organoids vs Tumor
surfme.3D_vs_Tumor <- get_CategoryExpressionPlot(results.3D_vs_Tumor$df, 
                                                 genes = surfme_diff$`3D`)
ggsave(file.path(date, plots_dir, "3D_vs_Tumor_surfme.png"),
       surfme.3D_vs_Tumor, "png", width = 10, height = 8, units = 'in')

