# ---------------------------------------------------------------------------- #
# -         Differential gene expression - Supplementary figures             - #
# ---------------------------------------------------------------------------- #
pca <- list()
# --- PCA plot: Patient 1 --- #
(pca[["593"]] <- make_pca(deseq.patients$`593`$rld, group = "condition", 
                         labs = c("Tumor tissue",# = "physioxia.Tumour",
                                  "3D spheroid (physioxia)", # = "physioxia.3D",
                                  "Cell culture (2D normoxia)", # = "normoxia.2D",
                                  "Cell culture (2D hypoxia)"), # = "hypoxia.2D"), 
                         cols = c("normoxia.2D" = "steelblue",
                                  "hypoxia.2D" = "salmon",
                                  "physioxia.3D" = "red",
                                  "physioxia.Tumour" = "darkred")))
# Save the plots
## as .svg file for publication
ggsave(file.path(date, plots_dir,"PCA_Patient1.svg"), plot = pca[["593"]],
       width = 12, height = 12, units = 'in', device = "svg")
## as .png file for presentation
ggsave(file.path(date, plots_dir,"PCA_Patient1.png"), plot = pca[["593"]],
       width = 12, height = 12, units = 'in', device = "png")

# --- PCA plot: Patient 2 --- #
(pca[["673"]] <- make_pca(deseq.patients$`673`$rld, group = "condition", 
                         labs = c("Tumor tissue",# = "physioxia.Tumour",
                                  "3D spheroid (physioxia)", # = "physioxia.3D",
                                  "Cell culture (2D normoxia)", # = "normoxia.2D",
                                  "Cell culture (2D hypoxia)"), # = "hypoxia.2D"), 
                         cols = c("normoxia.2D" = "steelblue",
                                  "hypoxia.2D" = "salmon",
                                  "physioxia.3D" = "red",
                                  "physioxia.Tumour" = "darkred")))
# Save the plots
## as .svg file for publication
ggsave(file.path(date, plots_dir,"PCA_Patient2.svg"), plot = pca[["673"]],
       width = 12, height = 12, units = 'in', device = "svg")
## as .png file for presentation
ggsave(file.path(date, plots_dir,"PCA_Patient2.png"), plot = pca[["673"]],
       width = 12, height = 12, units = 'in', device = "png")

# ------- Vulcano plot: Patient 1 - using 2DN as baseline expression --------- #
## 2.1 2D hypoxia vs 2D normoxia
(results.593$`2DN_vs_2DH`$plot <- make_vulcanplot(
  total = results.593$`2DN_vs_2DH`$df, 
  sig = results.593$`2DN_vs_2DH`$sig_df))
## 2.2 3D organoids vs 2D normoxia
(results.593$`2DN_vs_3D`$plot <- make_vulcanplot(
  total = results.593$`2DN_vs_3D`$df, 
  sig = results.593$`2DN_vs_3D`$sig_df))
## 2.3 Tumor samples vs 2D normoxia
(results.593$`2DN_vs_Tumor`$plot <- make_vulcanplot(
  total = results.593$`2DN_vs_Tumor`$df, 
  sig = results.593$`2DN_vs_Tumor`$sig_df))

## Merge the results with cowplot
(vulcanplots.593 <- merge_vulcanplots(
  p1 = results.593$`2DN_vs_2DH`$plot, # 2D hypoxia vs 2D normoxia
  p2 = results.593$`2DN_vs_3D`$plot, # 3D organoids vs 2D normoxia
  p3 = results.593$`2DN_vs_Tumor`$plot, # Tumor samples vs 2D normoxia
  titles = c("2D hypoxia vs 2D normoxia", 
             "3D organoids vs 2D normoxia", 
             "Tumor samples vs 2D normoxia")
))

# Save the plot
## as .svg file for publication
ggsave(file.path(date, plots_dir,"Patient593_total_vulcanplot.svg"), 
       plot = vulcanplots.593, width = 20, height = 8, units = 'in', device = "svg")
## as .png file for presentation
ggsave(file.path(date, plots_dir,"Patient593_total_vulcanplot.png"), 
       plot = vulcanplots.593, width = 20, height = 8, units = 'in', device = "png")

# ------- Vulcano plot: Patient 2 - using 2DN as baseline expression --------- #
## 2.1 2D hypoxia vs 2D normoxia
(results.673$`2DN_vs_2DH`$plot <- make_vulcanplot(
  total = results.673$`2DN_vs_2DH`$df, 
  sig = results.673$`2DN_vs_2DH`$sig_df))
## 2.2 3D organoids vs 2D normoxia
(results.673$`2DN_vs_3D`$plot <- make_vulcanplot(
  total = results.673$`2DN_vs_3D`$df, 
  sig = results.673$`2DN_vs_3D`$sig_df))
## 2.3 Tumor samples vs 2D normoxia
(results.673$`2DN_vs_Tumor`$plot <- make_vulcanplot(
  total = results.673$`2DN_vs_Tumor`$df, 
  sig = results.673$`2DN_vs_Tumor`$sig_df))

## Merge the results with cowplot
(vulcanplots.673 <- merge_vulcanplots(
  p1 = results.673$`2DN_vs_2DH`$plot, # 2D hypoxia vs 2D normoxia
  p2 = results.673$`2DN_vs_3D`$plot, # 3D organoids vs 2D normoxia
  p3 = results.673$`2DN_vs_Tumor`$plot, # Tumor samples vs 2D normoxia
  titles = c("2D hypoxia vs 2D normoxia", 
             "3D organoids vs 2D normoxia", 
             "Tumor samples vs 2D normoxia")
))

# Save the plot
## as .svg file for publication
ggsave(file.path(date, plots_dir,"Patient673_total_vulcanplot.svg"), 
       plot = vulcanplots.673, width = 20, height = 8, units = 'in', device = "svg")
## as .png file for presentation
ggsave(file.path(date, plots_dir,"Patient673_total_vulcanplot.png"), 
       plot = vulcanplots.673, width = 20, height = 8, units = 'in', device = "png")

# ------------------------------------------------------------------------------
# -         Figure 1.) COMPARING PATIENT 593 AND PATIENT 673                   -
# ------------------------------------------------------------------------------
### Venn diagram
LFC_range <- range(min(c(upset.593$log2FoldChange, upset.673$log2FoldChange)),
                   max(c(upset.593$log2FoldChange, upset.673$log2FoldChange)))

# ------- Venn diagram: Patient 1 - compare DEGS 2DH, 3D, Tumor vs 2DN ------- #
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
  + geom_venn_region(upset.593, sets = c("2DH","3D", "Tumor"), 
                     alpha = 0.3, show.legend = F) # add the regions
  + geom_venn_circle(upset.593, sets = c("2DH","3D", "Tumor"),
                     size = .5) # add the circles
  + scale_fill_venn_mix(upset.593, sets = c("2DH","3D", "Tumor"),
                        guide='none', highlight='none',
                        inactive_color='NA') # highlight the 3-way intersection
  + geom_venn_label_set(upset.593, outwards_adjust = 2,
                        sets = c("2DH","3D", "Tumor"),
                        fill = alpha("black", .35), 
                        aes(label = c("Tumor","3D", "2DH"))) 
  # add the set labels
  + geom_venn_label_region(
    upset.593, sets = c("2DH","3D", "Tumor"),
    aes(label=size), outwards_adjust=1.25, position=position_nudge(y=0.2)) 
  # add the number of genes in each region
  + theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )))

# save the plot
## as .svg file for publication
ggsave(file.path(date, plots_dir,"Patient593_venn_plot.svg"), venn_plot.593,
       "svg", width = 12, height = 6, units = 'in')
## as .png file for presentation
ggsave(file.path(date, plots_dir,"Patient593_venn_plot.png"), venn_plot.593,
       "png", width = 12, height = 6, units = 'in')


# ------- Venn diagram: Patient 2 - compare DEGS 2DH, 3D, Tumor vs 2DN ------- #
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
  + geom_venn_region(upset.673, sets = c("2DH","3D","Tumor"), 
                     alpha = 0.3, show.legend = F) # add the regions
  + geom_venn_circle(upset.673, sets = c("2DH","3D","Tumor"),
                     size = .5) # add the circles
  + scale_fill_venn_mix(upset.673, sets = c("2DH","3D","Tumor"),
                        guide='none', highlight = 'none',
                        inactive_color='NA') # highlight the 3-way intersection
  + geom_venn_label_set(upset.673, outwards_adjust = 2,
                        sets = c("2DH","3D","Tumor"),
                        fill = alpha("black", .35), 
                        aes(label = c("Tumor","3D","2DH")))) 
 # add the set labels
 + geom_venn_label_region(
   upset.673, sets = c("2DH","3D","Tumor"),
   aes(label=size), outwards_adjust=1.25, position=position_nudge(y=0.2)) 
 # add the number of genes in each region
 + theme(
   legend.title = element_text(size = 12),
   legend.text = element_text(size = 12)
 ))

# Save the plot
## as .svg file for publication
ggsave(file.path(date, plots_dir, "Patient673_venn_plot.svg"), venn_plot.673,
       "svg", width = 14, height = 10, units = 'in')
## as .png file for presentation
ggsave(file.path(date, plots_dir, "Patient673_venn_plot.png"), venn_plot.673,
       "png", width = 14, height = 10, units = 'in')

venn.legend <- get_legend(venn_plot.593)

venn.grid <- cowplot::plot_grid(venn_plot.593 + 
                     theme(legend.position = "none"),
                   venn_plot.673 + 
                     theme(legend.position = "none"),
                   venn.legend, 
                   labels = c("Patient #1","Patient #2"),
                   ncol = 3, rel_widths = c(1,1,0.3))

# Save the plot
## as .svg file for publication
ggsave(file.path(date, plots_dir, "complete_venn_plot.svg"), venn.grid,
       "svg", width = 20, height = 8, units = 'in')
## as .png file for presentation
ggsave(file.path(date, plots_dir, "complete_venn_plot.png"), venn.grid,
       "png", width = 20, height = 8, units = 'in')


# ------------------------------------------------------------------------------
# -         Figure 1.) ORA compare                             -
# ------------------------------------------------------------------------------
# 1.) Convert data frames to long format
(ORA.compare.plot <- compare_GO(ORA.593.compare, ORA.673.compare))

# Save the plot
## as .svg file for publication
ggsave(file.path(date, plots_dir,"ORA_total_compare.svg"), ORA.compare.plot,
       "svg", width = 20, height = 10, units = 'in')
## as .png file for presentation
ggsave(file.path(date, plots_dir,"ORA_total_compare.png"), ORA.compare.plot,
       "png", width = 20, height = 10, units = 'in')

compare_GO(GO.593.compare, GO.673.compare)
# ------------------------------------------------------------------------------
# -         Figure 1.) SURFME genes visualization                              -
# ------------------------------------------------------------------------------
# ------- SURFME filter: Patient 1 ------- # total: 1306
table(results.593$`2DN_vs_2DH`$df[which(results.593$`2DN_vs_2DH`$df$geneID %in% surfme$geneID),"significance"])
# UP: 15, DOWN 3, NO CHANGE: 1288, NOT EXPRESSED: 1831
table(results.593$`2DN_vs_3D`$df[which(results.593$`2DN_vs_3D`$df$geneID %in% surfme$geneID),"significance"])
# UP: 36, DOWN 25, NO CHANGE: 1245, NOT EXPRESSED: 1831
table(results.593$`2DN_vs_Tumor`$df[which(results.593$`2DN_vs_Tumor`$df$geneID %in% surfme$geneID),"significance"])
# UP: 54, DOWN 56, NO CHANGE: 1196, NOT EXPRESSED: 1831

# ------- SURFME filter: Patient 2 ------- # total: 1419
table(results.673$`2DN_vs_2DH`$df[which(results.673$`2DN_vs_2DH`$df$geneID %in% surfme$geneID),"significance"])
# UP: 13, DOWN 1, NO CHANGE: 1405, NOT EXPRESSED: 1718
table(results.673$`2DN_vs_3D`$df[which(results.673$`2DN_vs_3D`$df$geneID %in% surfme$geneID),"significance"])
# UP: 63, DOWN 97, NO CHANGE: 1259, NOT EXPRESSED: 1718
table(results.673$`2DN_vs_Tumor`$df[which(results.673$`2DN_vs_Tumor`$df$geneID %in% surfme$geneID),"significance"])
# UP: 78, DOWN 102, NO CHANGE: 1239, NOT EXPRESSED: 1718

SURFME.patients <- data.frame(
  "Trend" = c("UP","DOWN","NO CHANGE","NOT EXPRESSED"),
  "Patient1.2DN_vs_2DH" = c(15,3,1288,1831),
  "Patient1.2DN_vs_3D" = c(36,25,1245,1831),
  "Patient1.2DN_vs_Tumor" = c(54,56,1196,1831),
  "Patient2.2DN_vs_2DH" = c(13,1,1405,1718),
  "Patient2.2DN_vs_3D" = c(63,97,1259,1718),
  "Patient2.2DN_vs_Tumor" = c(78,102,1239,1718)
) %>% 
  pivot_longer(cols = -Trend, names_pattern = "(.*)\\.(.*)$", names_to = c("patient","condition"), values_to = "Count") %>% 
  dplyr::mutate(across(where(is.character), as.factor)) %>% 
  dplyr::mutate(Trend = factor(Trend, levels = c("NOT EXPRESSED","NO CHANGE","UP","DOWN")))

MF_colors <- c(
  # Binding Activities
  "binding" = "#00008A", "calcium ion binding" = "#0000FF", "ATP binding" = "#3399FF",
  "protein binding" = "#1591EA", "glycosaminoglycan binding" = "#00CCCC", 
  "MHC class II protein complex binding" = "#66FFFF",
  # Enzymatic and Catalytic Activities
  "protein kinase activity" = "#009900", "endopeptidase activity" = "#66CC00",
  # Receptor and Signaling Activities
  "signaling receptor binding" = "purple3", "G protein-coupled receptor activity" = "#990099",
  "transmembrane signaling receptor activity" = "#FF0077", 
  #Transport and Channel Activities
  "transmembrane transporter activity" = "#CC0000", 
  "monoatomic ion transmembrane transporter activity" = "#FF8000",
  #Other activities
  "other" = "#EBEBEB"
)

# ------- Patient 1 -------
(set.593 <- ggplot(data = subset(SURFME.patients, patient == "Patient1"), 
               aes(x = condition, y = Count, fill = Trend)) +
  #stacked bar plot
  geom_bar(width = .8, color = "black", stat = "identity", position = "stack") +
  # flip to horizontal
  coord_flip() + 
  scale_fill_manual(name = "Regulation trend",
                    values = c("UP" = "red", "DOWN" = "steelblue",
                               "NO CHANGE" = "lightgrey", "NOT EXPRESSED" = "white"),
                    labels = c("UP" = "Up-regulated", "DOWN" = "Down-regulated",
                               "NO CHANGE" = "No change", "NOT EXPRESSED" = "Not expressed")) +
  theme_minimal() +
  geom_label_repel(aes(label=paste0(Count," (",round((Count/3137)*100,2),"%)"),
                       fill = Trend), show.legend = F, size = 5,
                   position = position_stack(vjust = .5)) +
   labs(x = "", y = "# of genes") +
   theme(axis.text.y = element_text(size = 16),
         axis.text.x = element_text(size = 16),
         axis.title = element_text(size = 16, face = "bold"),
         legend.title = element_text(size = 16),
         legend.text = element_text(size = 14),
         legend.position = "right"))

(upset_MF_plot.593 <- upset(data = SURFME.593$upset_base,
                            intersect = c("2DN_vs_2DH","2DN_vs_3D","2DN_vs_Tumor"), 
                            name='Conditions',
                            mode='distinct',
                            min_degree = 1,
                            max_degree = 3,
                            min_size = 1,
                            base_annotations = list(
                              'Intersection size'=(
                                intersection_size(
                                  counts = T, size = 1,
                                  bar_number_threshold = 1,
                                  color = "#222222",
                                  mapping=aes(
                                    fill=TERM.parent))
                                + scale_fill_manual(values = MF_colors)
                                + guides(fill = guide_legend(ncol = 1))
                                + theme(legend.position = "right",
                                        legend.spacing.y = unit(10,'mm'),
                                        legend.direction = "vertical",
                                        axis.text = element_text(size=14),
                                        axis.title = element_text(size=16),
                                        legend.key.size = unit(10,'mm'),
                                        legend.title = element_text(size=16, face = "bold"),
                                        legend.text = element_text(size=14),
                                        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
                                + labs(y='# of genes', fill = "Molecular functions")
                              )
                            ),
                            themes=upset_modify_themes(
                              list(
                                'intersections_matrix'=theme(text=element_text(size=16)),
                                'overall_sizes'=theme(text=element_text(size=16))
                              )),
                            stripes=c('white','white','white')))

# ------- Patient 2 -------
(set.673 <- ggplot(data = subset(SURFME.patients, patient == "Patient2"), 
                   aes(x = condition, y = Count, fill = Trend)) +
   #stacked bar plot
   geom_bar(width = .8, color = "black", stat = "identity", position = "stack") +
   # flip to horizontal
   coord_flip() + 
   scale_fill_manual(name = "Regulation trend",
                     values = c("UP" = "red", "DOWN" = "steelblue",
                                "NO CHANGE" = "lightgrey", "NOT EXPRESSED" = "white"),
                     labels = c("UP" = "Up-regulated", "DOWN" = "Down-regulated",
                                "NO CHANGE" = "No change", "NOT EXPRESSED" = "Not expressed")) +
   theme_minimal() +
   geom_label_repel(aes(label=paste0(Count," (",round((Count/3137)*100,2),"%)"),
                        fill = Trend), show.legend = F, size = 5,
                    position = position_stack(vjust = .5))+
   labs(x = "", y = "# of genes") +
   theme(axis.text.y = element_text(size = 16),
         axis.text.x = element_text(size = 16),
         axis.title = element_text(size = 16, face = "bold"),
         legend.title = element_text(size = 16),
         legend.text = element_text(size = 14),
         legend.position = "right"))

(upset_MF_plot.673 <- upset(data = SURFME.673$upset_base,
                            intersect = c("2DN_vs_2DH","2DN_vs_3D","2DN_vs_Tumor"), 
                            name='Conditions',
                            mode='distinct',
                            min_degree = 1,
                            max_degree = 3,
                            min_size = 1,
                            base_annotations = list(
                              'Intersection size'=(
                                intersection_size(
                                  counts = T, size = 1,
                                  bar_number_threshold = 1,
                                  color = "#222222",
                                  mapping=aes(
                                    fill=TERM.parent))
                                + scale_fill_manual(values = MF_colors)
                                + guides(fill = guide_legend(ncol = 1))
                                + theme(legend.position = "right",
                                        legend.spacing.y = unit(10,'mm'),
                                        legend.direction = "vertical",
                                        axis.text = element_text(size=14),
                                        axis.title = element_text(size=16),
                                        legend.key.size = unit(10,'mm'),
                                        legend.title = element_text(size=16, face = "bold"),
                                        legend.text = element_text(size=14),
                                        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
                                + labs(y='# of genes', fill = "Molecular functions")
                              )
                            ),
                            themes=upset_modify_themes(
                              list(
                                'intersections_matrix'=theme(text=element_text(size=16)),
                                'overall_sizes'=theme(text=element_text(size=16))
                              )),
                            stripes=c('white','white','white')))

# Combine barplot and upset plot
#surfme.legend <- get_legend(upset_MF_plot.673)
surfme.grid <- cowplot::plot_grid(upset_MF_plot.593, upset_MF_plot.673,
                                  ncol = 2, rel_widths = c(1, 1))

ggsave(file.path(date, plots_dir, "surfme_MF_upset.png"), "png", 
       plot = surfme.grid, width = 20, height = 8, dpi = 300)

ggsave(file.path(date, plots_dir, "surfme_MF_upset.svg"), "svg", 
       plot = surfme.grid, width = 20, height = 8, units = 'in')

