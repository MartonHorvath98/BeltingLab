################################################################################
# 1.) Comparing the Uppsala cells and Lipid loaded cells   #
################################################################################

# GSEA on the GO terms from MSigDB
TOTAL.NES <- list(rbind(HGCC.GSEA$U3017$sig_df[,c("ID","Name", "NES")],
                        HGCC.GO$U3017$sig_df[,c("ID","Name", "NES")]),
                  rbind(HGCC.GSEA$U3047$sig_df[,c("ID","Name", "NES")],
                        HGCC.GO$U3047$sig_df[,c("ID","Name", "NES")]),
                  rbind(HGCC.GSEA$U3054$sig_df[,c("ID","Name", "NES")],
                        HGCC.GO$U3054$sig_df[,c("ID","Name", "NES")]),
                  rbind(CCLD.GSEA$sig_df[,c("ID","Name", "NES")],
                        CCLD.GO$sig_df[,c("ID","Name", "NES")]))
# Merge the data frames
TOTAL.NES <- merge.rec(TOTAL.NES, by = c("ID","Name"), 
                        suffixes = c("",""), all = T) %>% 
  setNames(., c("ID","Name", c("U3017","U3047","U3054","CCLD")))

TOTAL.FDR <- list(rbind(HGCC.GSEA$U3017$sig_df[,c("ID","Name", "p.adjust")],
                        HGCC.GO$U3017$sig_df[,c("ID","Name", "p.adjust")]),
                  rbind(HGCC.GSEA$U3047$sig_df[,c("ID","Name", "p.adjust")],
                        HGCC.GO$U3047$sig_df[,c("ID","Name", "p.adjust")]),
                  rbind(HGCC.GSEA$U3054$sig_df[,c("ID","Name", "p.adjust")],
                        HGCC.GO$U3054$sig_df[,c("ID","Name", "p.adjust")]),
                  rbind(CCLD.GSEA$sig_df[,c("ID","Name", "p.adjust")],
                        CCLD.GO$sig_df[,c("ID","Name", "p.adjust")]))
# Merge the data frames
TOTAL.FDR <- merge.rec(TOTAL.FDR, by = c("ID","Name"), 
                       suffixes = c("",""), all = T) %>% 
  setNames(., c("ID","Name", c("U3017.FDR","U3047.FDR","U3054.FDR","CCLD.FDR")))

# Extract the allocation of each gene
regions <- TOTAL.NES %>%
  # Pivot the data frame to long format along the expression values
  tidyr::pivot_longer(cols = !c(ID, Name), names_to = "Region", values_to = "Value") %>%
  # Filter out missing values
  dplyr::filter(!is.na(Value)) %>%
  # Group by geneID and concatenate the regions for genes in intersecting sets
  dplyr::group_by(ID, Name) %>%
  dplyr::summarise(Regions = paste(Region, collapse = "-")) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(Regions = as.factor(Regions))

# Merge the regions with the table
TOTAL.sets <- dplyr::inner_join(TOTAL.NES, TOTAL.FDR, by = c("Name", "ID")) %>%
  dplyr::inner_join(., regions, by = c("Name", "ID"))

# Calculate the number of genes in each subset
regions <- regions %>% 
  dplyr::group_by(Regions) %>% 
  dplyr::summarise(Size = n()) %>% 
  dplyr::pull(Size, name = Regions)


# Create a data frame for plotting the ellipses on the Venn diagram
region.matrix <- data.frame(
  labels = c("U3017","U3047","U3054","CCLD"),
  # Define the coordinates and dimensions of the ellipses: A, B, C, D
  x0=c(0.35, 0.5, 0.5, 0.65),
  y0=c(0.47, 0.57, 0.57, 0.47),
  a=c(0.35, 0.35, 0.35, 0.35),
  b=c(0.2, 0.15, 0.15, 0.2),
  angle=c(-10, -10, 10, 10),
  stringsAsFactors = T)
  
# Create a data frame for plotting the labels on the Venn diagram
label.matrix <- matrix(byrow = T, nrow = 15, ncol = 3,
                       # The labels equal the number of genes in each region
                       c(regions[['U3047']], 0.35, 0.75,
                         regions[['U3047-U3054']], 0.5, 0.69,
                         regions[['U3054']], 0.65, 0.75,
                         regions[['U3017-U3047']], 0.27, 0.65,
                         regions[['U3017-U3047-U3054']], 0.4, 0.6,
                         regions[['U3017-U3047-U3054-CCLD']], 0.5, 0.5,
                         regions[['U3047-U3054-CCLD']], 0.6, 0.6,
                         regions[['U3054-CCLD']], 0.73, 0.65,
                         regions[['U3017']], 0.15, 0.6,
                         regions[['U3017-U3054']], 0.3, 0.45,
                         regions[['U3017-U3054-CCLD']], 0.4, 0.4,
                         regions[['U3017-U3047-CCLD']], 0.6, 0.4,
                         regions[['U3047-CCLD']], 0.7, 0.45,
                         regions[['CCLD']], 0.85, 0.6,
                         regions[['U3017-CCLD']], 0.5, 0.28))
colnames(label.matrix) <- c("label", "x", "y")
row.names(label.matrix) <- names(regions)[c(10,12,14,4,6,7,13,15,2,8,9,5,11,1,3)]
label.matrix <- label.matrix %>%
  as.data.frame() %>% 
  tibble::rownames_to_column("Regions")

set.seed(0)
# Create a data frame for plotting the genes on the Venn diagram
point.matrix <- dplyr::inner_join(TOTAL.sets, label.matrix, 
                      by = c("Regions"))
point.matrix <- point.matrix %>% 
  dplyr::rowwise(.) %>% 
  # Add jitter to the coordinates for better visualization
  dplyr::mutate(
    # Determine the trend of gene expression
    trend = dplyr::case_when(
      all(dplyr::c_across(U3017:CCLD) > 0, na.rm = T) ~ 'UP',
      all(dplyr::c_across(U3017:CCLD) < 0, na.rm = T) ~ 'DOWN',
      TRUE ~ 'CHANGE'
    ),
    trend = as.factor(trend)
  ) %>% 
  dplyr::ungroup(.)

TOTAL.sets <- dplyr::inner_join(TOTAL.sets, point.matrix[,c("ID","trend")], 
                                by = c("ID"))

trend.matrix <- point.matrix %>% 
  dplyr::group_by(Regions, trend) %>% 
  dplyr::summarise(Size = n()) %>% 
  dplyr::inner_join(label.matrix[,-2], by = "Regions") %>% 
  tidyr::pivot_wider(names_from = trend, values_from = Size) %>% 
  dplyr::mutate_all(~ifelse(is.na(.),0,.))

  
# Create the Venn diagram
(TOTAL.venn <- ggplot() +
  # Clear the background and axes
  theme_dendro() +
  # Plot the ellipses
  ggforce::geom_ellipse(
    data = region.matrix, # ellipse coordinates
    aes(x0=x0,y0=y0,a=a,b=b, angle=angle, color = labels),
    fill = "white", linewidth = 1.2, alpha = .3, show.legend = F) +
  # Paint the ellipses
  scale_color_manual(
    values = c("U3017" = "pink",
               "U3047" = "salmon",
               "U3054" = "darkred",
               "CCLD" = "steelblue"),
    labels = c("U3017" = "U3017",
               "U3047" = "U3047",
               "U3054" = "U3054",
               "CCLD" = "LD vs. noLD")) +
  # Plot the number of genes in each region
  geom_scatterpie(data = trend.matrix, 
                  cols = c("UP","DOWN","CHANGE"), color = "grey25",
                  aes(x=x, y=y, group=Regions), pie_scale = 3) +
  
  # Plot the number of genes in each region
  geom_label(
    data = label.matrix, # Label coordinates
    aes(x = x, y = y, label = label),
    size = 6, alpha = 0.5, fill = 'white', color = 'grey15') + 
  # Paint the genes
  scale_fill_manual(
    name = "Trend",
    values = c("UP" = "red",
               "DOWN" = "blue",
               "CHANGE" = "orange"),
    labels = c("UP"="Activated", 
               "DOWN"="Inhibited",
               "CHANGE"="Opposing regulation")) +
    # Increase font size and adjust legend placement
    theme(text = element_text(size = 16),
          legend.key.spacing = unit(.5, 'cm'),
          legend.spacing = unit(2, 'cm')) + 
    # Expand the plot limits
    expand_limits(x = 0, y = .9) + 
    # Add region labels
    annotate('label', x = 0.1, y = 0.75, size = 8, label = 'U3017',
             color = "pink") +
    annotate('label', x = 0.25, y = 0.85, size = 8, label = 'U3047',
             color = "salmon") +
    annotate('label', x = 0.75, y = 0.85, size = 8, label = 'U3057',
             color = "darkred") + 
    annotate('label', x = 0.9, y = 0.75, size = 8, label = 'LD vs. noLD',
             color = "steelblue") +
    # Specify legend position
    theme(legend.position = 'bottom'))

# Save the plot
ggsave(file.path(wd, results_dir,"..",paste0("TOTAL_GO+pathway_venn_",Sys.Date(), ".png")),
       bg = "white", plot = TOTAL.venn, width = 12, height = 8)

TOTAL.sets <- TOTAL.sets %>% 
  dplyr::select(ID, Name, Regions, trend, U3017, U3017.FDR, U3047, U3047.FDR,
                U3054, U3054.FDR, CCLD, CCLD.FDR) %>%
  dplyr::rename(Trend = trend,
                U3017.NES = U3017, U3047.NES = U3047,
                U3054.NES = U3054, CCLD.NES = CCLD) %>%
  dplyr::group_split(Regions, .keep = T)
names(TOTAL.sets) <- names(regions)

# Save the table
sapply(names(TOTAL.sets), function(x){
  openxlsx::write.xlsx(TOTAL.sets[[x]], 
                       file.path(wd, results_dir, "..", 
                                 paste0("TOTAL_GO+pathway_venn_", x, "_", date,".xlsx")))
})

rm(list = c(ls(pattern = "matrix"), "regions", "TOTAL.NES", "TOTAL.FDR"))

################################################################################
# 3. Visualize selected enrichment scores                                      #
################################################################################

library(enrichplot)
palette = c("EXTRACELLULAR_MATRIX_ORGANIZATION" = "#380186",
            "GLYCOSAMINOGLYCAN_METABOLISM" = "#8D00FA",
            "PROTEOGLYCAN_METABOLIC_PROCESS" = "#E70041",
            "CHONDROITIN_SULFATE_DERMATAN_SULFATE_METABOLISM" = "#FF08FF",
            "FAT_CELL_DIFFERENTIATION" = "#000000")
pathway.ranks = which(CCLD.GSEA$df$Name %in% names(palette))
go.ranks = which(CCLD.GO$df$Name %in% names(palette))

CCLD.enrichplot <- list()
CCLD.enrichplot$runningScores  <- do.call(rbind, c(lapply(pathway.ranks, gsInfo, object = CCLD.GSEA$gsea),
                                                   lapply(go.ranks, gsInfo, object = CCLD.GO$gsea)))
CCLD.enrichplot$runningScores$Description <- factor(
  gsub(x = CCLD.enrichplot$runningScores$Description, pattern = "REACTOME_|GOBP_", replacement =  ""),
  levels = names(palette))

p1 <- plotRunningScore(.df = CCLD.enrichplot$runningScores, 
                       .x = "x", .y = "runningScore", 
                       .color = "Description", .palette = palette)

p2 <- plotGeneRank(.df = CCLD.enrichplot$runningScores, 
             .x = "x", .facet = "Description~.",
             .color = "Description", .palette = palette)


CCLD.enrichplot$table <- rbind(CCLD.GSEA$df[pathway.ranks,],
                 CCLD.GO$df[go.ranks,]) %>% 
  dplyr::arrange(order(match(Name, names(palette)))) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("Name") %>%
  dplyr::mutate(
    NES = round(NES, 4),
    FDR = round(p.adjust, 4)) %>%
  dplyr::mutate(FDR = case_when(
    FDR < 0.001 ~ paste("<0.001", "(***)"),
    FDR < 0.01 ~ paste(as.character(FDR), "(**)"),
    FDR < 0.05 ~ paste(as.character(FDR), "(*)"),
    TRUE ~ as.character(FDR),
  )) %>% 
  dplyr::select(c(5,13))

table_grob <- tableGrob(
  CCLD.enrichplot$table, 
  theme = ttheme_minimal(
    base_size = 14,
    core = list(
      fg_params = list(hjust = 0.5, x = 0.5, col = palette)),
    rowhead = list(
      fg_params = list(hjust = 0, x = 0,col = palette[c(5,1,2,3,4)]))
    )) 

CCLD.enrichplot$plot <- cowplot::plot_grid(
  p1, "", p2, table_grob,
  byrow = T, nrow = 2, ncol = 2, 
  rel_heights = c(1.5,1), rel_widths =  c(1,2),
  margins = c(0.5, 0.5, 0.5, 0.5))

ggsave(file.path(results_dir, date, "plots", "CCLD_GSEA.svg"), bg = "white", device = "svg",
       plot = CCLD.enrichplot$plot, width = 11.7, height = 8.3, units = "in")

interest_genes <- read.csv("genes-of-interest.txt", header = T,
                           sep = "\t", stringsAsFactors = T)

heatmap_data <- interest_genes %>% 
  dplyr::inner_join(., CC_LDvsnoLD[,3:5], by = c("Gene.name" = "SYMBOL")) %>%
  dplyr::inner_join(., HGCC.df$U3017[,c(1,2)], by = c("Gene.name" = "symbol")) %>% 
  dplyr::inner_join(., HGCC.df$U3047[,c(1,2)], by = c("Gene.name" = "symbol")) %>%
  dplyr::inner_join(., HGCC.df$U3054[,c(1,2)], by = c("Gene.name" = "symbol")) %>% 
  dplyr::rename_at(vars(contains("log2")), ~c("LDvsnoLD", "U3017", "U3047", "U3054")) %>%
  dplyr::mutate(Category = factor(Category, 
                                      levels = c("Acidosis/Hypoxia", "CSPG core", "CSPG Biosynthesis- GAG linker",
                                              "CSPG Biosynthesis","Lipid droplet"),
                                      labels = c("Acidosis/Hypoxia", "CSPG core", "GAG linker",
                                              "CSPG Biosynthesis","Lipid droplet"))) %>% 
  # dplyr::mutate(avgFC = rowMeans(.[,c(4:7)])) %>% 
  # dplyr::group_by(Category) %>%
  # dplyr::mutate(rank = rank(avgFC, ties.method = "min")) %>% 
  # dplyr::arrange(Category, desc(rank)) %>% 
  # dplyr::ungroup() %>% 
  dplyr::select(Gene.name, Category, U3017, U3047, U3054, LDvsnoLD)

row_order <- match(c("CA9","VEGFA","CA12","NCAN","CSPG5","VCAN","BCAN","BGN",
                     "ACAN","DCN","CSPG4","CD44","XYLT1","XYLT2","CSGALNACT1",
                     "CHSY1","CHSY3","CHPF","CHPF2","DSEL","HILPDA", "FABP4",
                     "G0S2","PLIN2","PLIN1","FABP5"),
                   heatmap_data$Gene.name)
col_order <- order(c(as.numeric(proxy::dist(rbind(heatmap_matrix[,1, drop=T],
                                                  heatmap_matrix[,4, drop = T]),
                                            method = "euclidean", pairwise = T)),
                     as.numeric(proxy::dist(rbind(heatmap_matrix[,2, drop=T],
                                                  heatmap_matrix[,4, drop = T]),
                                            method = "euclidean", pairwise = T)),
                     as.numeric(proxy::dist(rbind(heatmap_matrix[,3, drop = T],
                                                  heatmap_matrix[,4, drop = T]),
                                            method = "euclidean", pairwise = T)),
                     as.numeric(proxy::dist(rbind(heatmap_matrix[,4, drop =T],
                                                  heatmap_matrix[,4, drop=T]),
                                            method = "euclidean", pairwise = T))))

heatmap_data <- heatmap_data[row_order,]
heatmap_matrix <- heatmap_data[,c(3:6)]
heatmap_matrix <- as.matrix(heatmap_matrix[,col_order])
       
row.names(heatmap_matrix) <- heatmap_data$Gene.name

library(ComplexHeatmap)
library(circlize)
library(gridtext)

cat_annot = rowAnnotation(
  Category = heatmap_data$Category,
  col = list(Category = c("Acidosis/Hypoxia" = "darkgreen",
                          "CSPG core" = "turquoise",
                          "GAG linker" = "#8D00FA",
                          "CSPG Biosynthesis" = "#FF08FF",
                          "Lipid droplet" = "white")
            
  ),
  gp = gpar(col = "black", fontsize = 14, family = "arial"),
  show_annotation_name = F,
  annotation_legend_param = list(border = "black")
)

matrix_col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
(ht = Heatmap(heatmap_matrix, 
        name = "Log2 (Fold change)",
        col = matrix_col,
        
        cluster_rows = F,
        cluster_columns = F,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 14, family = "arial"),
        row_gap = unit(5, "mm"),
        column_names_rot = 0,
        column_names_centered = T,
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 14, family = "arial"),
        
        column_labels = c("LDvsnoLD" = "", "U3047" = "U3047",
                       "U3054" = "U3054", "U3017"="U3017"),
        column_split = factor(c("LD vs. noLD", rep("2D vs. 3D", 3)),
                           levels = c("LD vs. noLD", "2D vs. 3D")),
        
        right_annotation = cat_annot, 
        heatmap_legend_param = list(direction = "vertical"),
        border_gp = gpar(col = "black", lty = 2)))

svg(file.path("..", results_dir, "heatmap.svg"), width = 12, height = 8)#,
#    res = 300, units = "in")
draw(ht, merge_legend = T,
     heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

################################################################################
# 4. Vulcano visualization of the shared terms and pathways                    #
################################################################################
