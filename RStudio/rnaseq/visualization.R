##############################################
# 1.) Visualise DESeq2 results               #
##############################################
# 1. PCA plot
# ---------------- Total results -----------------------------------------------
(pca <- make_pca(deseq$rld, group = "condition", 
                 labs = c("Tumor tissue",# = "physioxia.Tumour",
                          "3D organoids (physioxia)", # = "physioxia.3D",
                          "Cell culture (2D normoxia)", # = "normoxia.2D",
                          "Cell culture (2D hypoxia)"), # = "hypoxia.2D"), 
                cols = c("normoxia.2D" = "steelblue",
                         "hypoxia.2D" = "salmon",
                         "physioxia.3D" = "red",
                         "physioxia.Tumour" = "darkred")))
# save the plot
ggsave(paste(file.path(date, plots_dir),"pca.png",sep="/"),
       plot = pca, width = 10, height = 11, units = 'in')

# ---------------- In-vitro results --------------------------------------------
(invitro.pca <- make_pca(invitro.deseq$rld, group = "condition", 
                    labs = c("Cell culture (2D normoxia)", # = "normoxia.2D",
                             "Cell culture (2D hypoxia)",
                             "3D organoids (physioxia)"), # = "hypoxia.2D"), 
                    cols = c("normoxia.2D" = "steelblue",
                             "hypoxia.2D" = "salmon",
                             "physioxia.3D" = "darkred")))


# 2. Vulcano plots
# ---------------- Total results -----------------------------------------------
## 2.1 Tumor samples vs 2D normoxia
(results.2DN_vs_Tumor$plot <- make_vulcanplot(
  total = results.2DN_vs_Tumor$df, 
  sig = results.2DN_vs_Tumor$sig_df))

# save the plot
ggsave(file.path(date, plots_dir, "Tumor_vs_2DN_vulcan.png"), 
       plot = results.2DN_vs_Tumor$plot, width = 10, height = 8, units = 'in')

## 2.2 Tumor samples vs 2D hypoxia
(results.2DH_vs_Tumor$plot <- make_vulcanplot(
  total = results.2DH_vs_Tumor$df, 
  sig = results.2DH_vs_Tumor$sig_df))

# save the plot
ggsave(file.path(date, plots_dir, "Tumor_vs_2DH_vulcan.png"), 
       plot = results.2DH_vs_Tumor$plot, width = 10, height = 8, units = 'in')

## 2.3 Tumor samples vs 3D organoids
(results.3D_vs_Tumor$plot <- make_vulcanplot(
  total = results.3D_vs_Tumor$df, 
  sig = results.3D_vs_Tumor$sig_df))

# save the plot
ggsave(file.path(date, plots_dir, "Tumor_vs_3D_vulcan.png"), 
       plot = results.3D_vs_Tumor$plot, width = 10, height = 8, units = 'in')

# ---------------- In-vitro results --------------------------------------------
## 2.4 2D hypoxia vs 2D normoxia
(results.2DH_vs_2DN$plot <- make_vulcanplot(
  total = results.2DH_vs_2DN$df, 
  sig = results.2DH_vs_2DN$sig_df))

# save the plot
ggsave(file.path(date, plots_dir, "2DH_vs_2DN_vulcan.png"), 
       plot = results.2DH_vs_2DN$plot, width = 10, height = 8, units = 'in')

## 2.5 3D organoids vs 2D normoxia
(results.3D_vs_2DN$plot <- make_vulcanplot(
  total = results.3D_vs_2DN$df, 
  sig = results.3D_vs_2DN$sig_df))

# save the plot
ggsave(file.path(date, plots_dir, "3D_vs_2DN_vulcan.png"), 
       plot = results.3D_vs_2DN$plot, width = 10, height = 8, units = 'in')

##############################################
# 2.) Visualise KEGG (ORA) results           #
##############################################
# ---------------- Total results -----------------------------------------------
## Create plots for visual analysis
(KEGG.2DN_vs_Tumor$plot <- make_dotplot(KEGG.2DN_vs_Tumor$df, 
                                        KEGG.2DN_vs_Tumor$df$Count, 
                                        type = "KEGG"))
ggsave(file.path(date, plots_dir,"2DN_vs_Tumor_KEGG.png"), 
       plot = KEGG.2DN_vs_Tumor$plot, width = 10, height = 8, units = 'in')

(KEGG.2DH_vs_Tumor$plot <- make_dotplot(KEGG.2DH_vs_Tumor$df, 
                                        KEGG.2DH_vs_Tumor$df$Count, 
                                        type = "KEGG"))
ggsave(file.path(date, plots_dir,"2DH_vs_Tumor_KEGG.png"),
       plot = KEGG.2DH_vs_Tumor$plot, width = 10, height = 8, units = 'in')

(KEGG.3D_vs_Tumor$plot <- make_dotplot(KEGG.3D_vs_Tumor$df, 
                                       KEGG.3D_vs_Tumor$df$Count, 
                                       type = "KEGG"))
ggsave(file.path(date, plots_dir,"3D_vs_Tumor_KEGG.png"),
       plot = KEGG.3D_vs_Tumor$plot, width = 10, height = 8, units = 'in')

# ---------------- In-vitro results --------------------------------------------
## Create plots for visual analysis
(KEGG.2DH_vs_2DN$plot <- make_dotplot(KEGG.2DH_vs_2DN$df, 
                                      KEGG.2DH_vs_2DN$df$Count, 
                                      type = "KEGG"))
ggsave(file.path(date, plots_dir,"2DH_vs_2DN_KEGG.png"),
       plot = KEGG.2DH_vs_2DN$plot, width = 10, height = 8, units = 'in')

(KEGG.3D_vs_2DN$plot <- make_dotplot(KEGG.3D_vs_2DN$df, 
                                     KEGG.3D_vs_2DN$df$Count, 
                                     type = "KEGG"))
ggsave(file.path(date, plots_dir,"3D_vs_2DN_KEGG.png"),
       plot = KEGG.3D_vs_2DN$plot, width = 10, height = 8, units = 'in')

##############################################
# 3.) Visualise GO (ORA) results             #
##############################################
# ---------------- Total results -----------------------------------------------
## Create plots for visual analysis
(GO.2DN_vs_Tumor$plot <- make_dotplot(GO.2DN_vs_Tumor$df, 
                                      GO.2DN_vs_Tumor$df$Count, 
                                      type = "GO"))
ggsave(file.path(date, plots_dir,"2DN_vs_Tumor_GO.png"),
       plot = GO.2DN_vs_Tumor$plot, width = 10, height = 8, units = 'in')

(GO.2DH_vs_Tumor$plot <- make_dotplot(GO.2DH_vs_Tumor$df, 
                                      GO.2DH_vs_Tumor$df$Count, 
                                      type = "GO"))
ggsave(file.path(date, plots_dir,"2DH_vs_Tumor_GO.png"),
       plot = GO.2DH_vs_Tumor$plot, width = 10, height = 8, units = 'in')

(GO.3D_vs_Tumor$plot <- make_dotplot(GO.3D_vs_Tumor$df, 
                                     GO.3D_vs_Tumor$df$Count, 
                                     type = "GO"))
ggsave(file.path(date, plots_dir,"3D_vs_Tumor_GO.png"),
       plot = GO.3D_vs_Tumor$plot, width = 10, height = 8, units = 'in')

# ---------------- In-vitro results --------------------------------------------
## Create plots for visual analysis
(GO.2DH_vs_2DN$plot <- make_dotplot(GO.2DH_vs_2DN$df, 
                                    GO.2DH_vs_2DN$df$Count, 
                                    type = "GO"))
ggsave(file.path(date, plots_dir,"2DH_vs_2DN_GO.png"),
       plot = GO.2DH_vs_2DN$plot, width = 10, height = 8, units = 'in')

(GO.3D_vs_2DN$plot <- make_dotplot(GO.3D_vs_2DN$df, 
                                   GO.3D_vs_2DN$df$Count, 
                                   type = "GO"))
ggsave(file.path(date, plots_dir,"3D_vs_2DN_GO.png"),
       plot = GO.3D_vs_2DN$plot, width = 10, height = 8, units = 'in')

##############################################
# 4.) GO Semantic similarity visual analysis #
##############################################
# ---------------- Total results -----------------------------------------------
## 1.1 2D normoxia vs Tumor - Biological Process
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

# ---------------- In-vitro results --------------------------------------------
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

##############################################
# 5.) Compound GO plots                      #
##############################################
# ---------------- Total results -----------------------------------------------
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

# ---------------- In-vitro results --------------------------------------------
## 4. 2D hypoxia vs 2D normoxia
(COMPPLOT.2DH_vs_2DN <- (compound_GOplot(COMP.2DH_vs_2DN)))
ggsave(file.path(date, plots_dir,"2DH_vs_2DN_compound_GO.png"),
       COMPPLOT.2DH_vs_2DN, "png", width = 20, height = 12, units = 'in')

## 5. 3D organoids vs 2D normoxia
(COMPPLOT.3D_vs_2DN <- (compound_GOplot(COMP.3D_vs_2DN)))
ggsave(file.path(date, plots_dir,"3D_vs_2DN_compound_GO.png"),
       COMPPLOT.3D_vs_2DN, "png", width = 20, height = 12, units = 'in')
################################################
# 6.) Venn diagram for DEGs                    #
################################################
# ---------------- Total results -----------------------------------------------
(venn_plot <- ( # make the plot
  ggplot(upset) 
  + coord_fixed()
  + theme_void() # remove the background
  + geom_point(aes(x = x, y = y, colour = log2FoldChange), size = 2) # add the genes
  + scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    na.value = "grey50", guide = "colourbar", aesthetics = "colour"
  ) # add a color gradient to the points
  + geom_venn_region(upset, sets = c("2DN", "2DH", "3D"), 
                     alpha = 0.3, show.legend = F) # add the regions
  + geom_venn_circle(upset, sets = c("2DN", "2DH", "3D"),
                     size = .5) # add the circles
  + scale_fill_venn_mix(upset, sets = c("2DN", "2DH", "3D"),
                        guide='none', highlight=c("2DN-2DH-3D"),
                        inactive_color='NA') # highlight the 3-way intersection
  + geom_venn_label_set(upset, outwards_adjust = 2,
                        sets = c("2DN", "2DH", "3D"),
                        fill = alpha("black", .35), 
                        aes(label = c("2DN", "2DH", "3D")))) 
  # add the set labels
  + geom_venn_label_region(
    upset, sets = c("2DN", "2DH", "3D"),
    aes(label=size), outwards_adjust=1.25, position=position_nudge(y=0.2)) 
  # add the number of genes in each region
  + theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  ))

# save the plot
ggsave(file.path(date, plots_dir, "venn_plot.png"), venn_plot,
       "png", width = 12, height = 6, units = 'in')
