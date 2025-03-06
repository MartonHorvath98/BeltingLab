################################################################################
# 1.) Comparing the Uppsala cells and Lipid loaded cells                       #
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
    annotate('label', x = 0.75, y = 0.85, size = 8, label = 'U3054',
             color = "darkred") + 
    annotate('label', x = 0.9, y = 0.75, size = 8, label = 'LD vs. noLD',
             color = "steelblue") +
    # Specify legend position
    theme(legend.position = 'bottom'))

# Save the plot
results_dir <- "Results/HGCC"
ggsave(file.path(results_dir,date,"plots",paste0("TOTAL_GO+pathway_venn_",Sys.Date(), ".png")),
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

(CCLD.enrichplot$plot <- cowplot::plot_grid(
  p1, p2, #table_grob,
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r", #rel_widths =  c(1,2),
  margins = c(0.5, 0.5, 0.5, 0.5)))

ggsave(file.path(results_dir, date, "plots", "CCLD_GSEA.svg"), bg = "white", device = "svg",
       plot = CCLD.enrichplot$plot, width = 6, height = 4.45, units = "in")

interest_genes <- read.csv("data/genes-of-interest.txt", header = T,
                           sep = "\t", stringsAsFactors = T)
interest_genes <- c("CSGALNACT1","CSGALNACT2","CHSY1","CHSY3","CHPF","CHPF2",
                    "DCN","BGN","DSE","DSEL","CA9","HIF1A","G0S2","C7orf68")

TOTAL.heatmap <- list()
TOTAL.heatmap$data <- interest_genes %>% 
  dplyr::inner_join(., CCLD.df[,1:2], by = c("Gene.name" = "Symbol")) %>%
  dplyr::inner_join(., HGCC.deg$U3017[,c(2,4)], by = c("Gene.name" = "Symbol")) %>% 
  dplyr::inner_join(., HGCC.deg$U3047[,c(2,4)], by = c("Gene.name" = "Symbol")) %>%
  dplyr::inner_join(., HGCC.deg$U3054[,c(2,4)], by = c("Gene.name" = "Symbol")) %>% 
  dplyr::rename_at(vars(contains("log2")), ~c("LDvsnoLD", "U3017", "U3047", "U3054")) %>%
  dplyr::mutate(Category = factor(Category, 
                                      levels = c("Acidosis/Hypoxia", "CSPG core", "CSPG Biosynthesis- GAG linker",
                                              "CSPG Biosynthesis","Lipid droplet"),
                                      labels = c("Acidosis/Hypoxia", "CSPG core", "GAG linker",
                                              "CSPG Biosynthesis","Lipid droplet"))) %>% 
  dplyr::select(Gene.name, Category, U3017, U3047, U3054, LDvsnoLD)

TOTAL.heatmap$matrix <- as.matrix(TOTAL.heatmap$data[,c(3:6)])

row_order <- match(c("CA9","VEGFA","CA12","NCAN","CSPG5","VCAN","BCAN","BGN",
                     "ACAN","DCN","CSPG4","CD44","XYLT1","XYLT2","CSGALNACT1",
                     "CHSY1","CHSY3","CHPF","CHPF2","DSEL","HILPDA", "FABP4",
                     "G0S2","PLIN2","PLIN1","FABP5"),
                   TOTAL.heatmap$data$Gene.name)
col_order <- order(c(
  as.numeric(proxy::dist(rbind(TOTAL.heatmap$matrix[,1, drop=T],
                               TOTAL.heatmap$matrix[,4, drop = T]),
                         method = "euclidean", pairwise = T)),
  as.numeric(proxy::dist(rbind(TOTAL.heatmap$matrix[,2, drop=T],
                               TOTAL.heatmap$matrix[,4, drop = T]),
                         method = "euclidean", pairwise = T)),
  as.numeric(proxy::dist(rbind(TOTAL.heatmap$matrix[,3, drop = T],
                               TOTAL.heatmap$matrix[,4, drop = T]),
                         method = "euclidean", pairwise = T)),
  as.numeric(proxy::dist(rbind(TOTAL.heatmap$matrix[,4, drop =T],
                               TOTAL.heatmap$matrix[,4, drop=T]),
                         method = "euclidean", pairwise = T))))


row.names(TOTAL.heatmap$matrix) <- TOTAL.heatmap$data$Gene.name
TOTAL.heatmap$matrix <- TOTAL.heatmap$matrix[row_order,col_order]

cat_annot = rowAnnotation(
  Category = TOTAL.heatmap$data$Category[row_order],
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
(TOTAL.heatmap$plot = Heatmap(TOTAL.heatmap$matrix, # data matrix
                              
                              ### row settings
                              cluster_rows = F,
                              row_gap = unit(5, "mm"),
                              row_names_side = "left",
                              row_names_gp = gpar(fontsize = 14,
                                                  family = "arial"),
                              
                              ### column settings
                              cluster_columns = F,
                              column_names_rot = 0,
                              column_names_centered = T,
                              column_names_side = "top",
                              column_names_gp = gpar(fontsize = 14, 
                                                     family = "arial"),
                              column_labels = c("LDvsnoLD" = "", 
                                                "U3047" = "U3047",
                                                "U3054" = "U3054", 
                                                "U3017"="U3017"),
                              column_split = factor(c("LD vs. noLD", rep("2D vs. 3D", 3)),
                                                    levels = c("LD vs. noLD", "2D vs. 3D")),
                              
                              ### legend settings
                              col = matrix_col, # color settings
                              name = "Log2 (Fold change)", # color legend title
                              right_annotation = cat_annot, 
                              heatmap_legend_param = list(direction = "vertical"),
                              border_gp = gpar(col = "black", lty = 2)))

svg(file.path(results_dir,"..", "category_heatmap.svg"), width = 12, height = 8)
draw(TOTAL.heatmap$plot , merge_legend = T,
     heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

################################################################################
# 4.) The comparison of AA vs CA vs Hypoxia in U87 cells                       #
################################################################################
# 1. Load the GSEA results 
# GSEA on the GO terms from MSigDB
U87.NES <- list(rbind(U87.GSEA$`sel_pH647-control_sel`$sig_df[,c("ID","Name", "NES")],
                      U87.GO$`sel_pH647-control_sel`$sig_df[,c("ID","Name", "NES")]),
                rbind(U87.GSEA$`acu_pH68-control_acu`$sig_df[,c("ID","Name", "NES")],
                      U87.GO$`acu_pH68-control_acu`$sig_df[,c("ID","Name", "NES")]),
                rbind(U87.GSEA$`acu_pH64-control_acu`$sig_df[,c("ID","Name", "NES")],
                      U87.GO$`acu_pH64-control_acu`$sig_df[,c("ID","Name", "NES")]),
                rbind(U87.GSEA$`hypoxia-control_nox`$sig_df[,c("ID","Name", "NES")],
                      U87.GO$`hypoxia-control_nox`$sig_df[,c("ID","Name", "NES")]))
# Merge the data frames
U87.NES <- merge.rec(U87.NES, by = c("ID","Name"), 
                       suffixes = c("",""), all = T) %>% 
  setNames(., c("ID","Name", c("pH6.4_CA","pH6.4_AA","pH6.8_AA","Hypoxia")))

U87.FDR <- list(rbind(U87.GSEA$`sel_pH647-control_sel`$sig_df[,c("ID","Name", "p.adjust")],
                     U87.GO$`sel_pH647-control_sel`$sig_df[,c("ID","Name", "p.adjust")]),
               rbind(U87.GSEA$`acu_pH68-control_acu`$sig_df[,c("ID","Name", "p.adjust")],
                     U87.GO$`acu_pH68-control_acu`$sig_df[,c("ID","Name", "p.adjust")]),
               rbind(U87.GSEA$`acu_pH64-control_acu`$sig_df[,c("ID","Name", "p.adjust")],
                     U87.GO$`acu_pH64-control_acu`$sig_df[,c("ID","Name", "p.adjust")]),
               rbind(U87.GSEA$`hypoxia-control_nox`$sig_df[,c("ID","Name", "p.adjust")],
                     U87.GO$`hypoxia-control_nox`$sig_df[,c("ID","Name", "p.adjust")]))
# Merge the data frames
U87.FDR <- merge.rec(U87.FDR, by = c("ID","Name"), 
                       suffixes = c("",""), all = T) %>% 
  setNames(., c("ID","Name", c("pH6.4_CA","pH6.4_AA","pH6.8_AA","Hypoxia")))

# Extract the allocation of each gene
regions <- U87.NES %>%
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
U87.sets <- dplyr::inner_join(U87.NES, U87.FDR, by = c("Name", "ID"), suffix = c(".NES",".FDR")) %>%
  dplyr::inner_join(., regions, by = c("Name", "ID"))

# Calculate the number of genes in each subset
regions <- regions %>% 
  dplyr::group_by(Regions) %>% 
  dplyr::summarise(Size = n()) %>% 
  dplyr::pull(Size, name = Regions)


# Create a data frame for plotting the ellipses on the Venn diagram
region.matrix <- data.frame(
  labels = c("pH6.4_CA","pH6.4_AA","pH6.8_AA","Hypoxia"),
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
                       c(regions[['pH6.4_AA']], 0.35, 0.75,
                         regions[['pH6.4_AA-pH6.8_AA']], 0.5, 0.69,
                         regions[['pH6.8_AA']], 0.65, 0.75,
                         regions[['pH6.4_CA-pH6.4_AA']], 0.27, 0.65,
                         regions[['pH6.4_CA-pH6.4_AA-pH6.8_AA']], 0.4, 0.6,
                         regions[['pH6.4_CA-pH6.4_AA-pH6.8_AA-Hypoxia']], 0.5, 0.5,
                         regions[['pH6.4_AA-pH6.8_AA-Hypoxia']], 0.6, 0.6,
                         regions[['pH6.8_AA-Hypoxia']], 0.73, 0.65,
                         regions[['pH6.4_CA']], 0.15, 0.6,
                         regions[['pH6.4_CA-pH6.8_AA']], 0.3, 0.45,
                         regions[['pH6.4_CA-pH6.8_AA-Hypoxia']], 0.4, 0.4,
                         regions[['pH6.4_CA-pH6.4_AA-Hypoxia']], 0.6, 0.4,
                         regions[['pH6.4_AA-Hypoxia']], 0.7, 0.45,
                         regions[['Hypoxia']], 0.85, 0.6,
                         regions[['pH6.4_CA-Hypoxia']], 0.5, 0.28))
colnames(label.matrix) <- c("label", "x", "y")
row.names(label.matrix) <- names(regions)[c(2,4,14,8,10,11,5,15,6,12,13,9,3,1,7)]
label.matrix <- label.matrix %>%
  as.data.frame() %>% 
  tibble::rownames_to_column("Regions")

# Create a data frame for plotting the genes on the Venn diagram
U87.sets <- U87.sets %>% 
  dplyr::rowwise(.) %>% 
  # Add jitter to the coordinates for better visualization
  dplyr::mutate(
    # Determine the trend of gene expression
    trend = dplyr::case_when(
      all(dplyr::c_across(`pH6.4_CA.NES`:`Hypoxia.NES`) > 0, na.rm = T) ~ 'UP',
      all(dplyr::c_across(`pH6.4_CA.NES`:`Hypoxia.NES`) < 0, na.rm = T) ~ 'DOWN',
      TRUE ~ 'CHANGE'
    ),
    trend = as.factor(trend)
  ) %>% 
  dplyr::ungroup(.)

trend.matrix <- U87.sets %>% 
  dplyr::group_by(Regions, trend) %>% 
  dplyr::summarise(Size = n()) %>% 
  dplyr::inner_join(label.matrix[,-2], by = "Regions") %>% 
  tidyr::pivot_wider(names_from = trend, values_from = Size) %>% 
  dplyr::mutate_all(~ifelse(is.na(.),0,.))


# Create the Venn diagram
(U87.venn <- ggplot() +
    # Clear the background and axes
    theme_dendro() +
    # Plot the ellipses
    ggforce::geom_ellipse(
      data = region.matrix, # ellipse coordinates
      aes(x0=x0,y0=y0,a=a,b=b, angle=angle, color = labels),
      fill = "white", linewidth = 1.2, alpha = .3, show.legend = F) +
    # Paint the ellipses
    scale_color_manual(
      values = c("pH6.8_AA" = "pink",
                 "pH6.4_AA" = "brown1",
                 "pH6.4_CA" = "darkred",
                 "Hypoxia" = "purple"),
      labels = c("pH6.8_AA" = "Acute acidosis (pH6.8)",
                 "pH6.4_AA" = "Acute acidosis (pH6.4)",
                 "pH6.4_CA" = "Chronic acidosis (pH6.4)",
                 "Hypoxia" = "Hypoxia")) +
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
    annotate('label', x = 0.75, y = 0.85, size = 8, label = 'Acute acidosis (pH6.8)',
             color = "pink") +
    annotate('label', x = 0.25, y = 0.85, size = 8, label = 'Acute acidosis (pH6.4)',
             color = "brown1") +
    annotate('label', x = 0.1, y = 0.75, size = 8, label = 'Chronic acidosis (pH6.4)',
             color = "darkred") + 
    annotate('label', x = 0.9, y = 0.75, size = 8, label = 'Hypoxia',
             color = "purple") +
    # Specify legend position
    theme(legend.position = 'bottom'))

# Save the plot
ggsave(file.path(wd, results_dir,date, "plots","U87_GO+pathway_venn.png"),
       bg = "white", plot = U87.venn, width = 12, height = 8)

U87.sets <- U87.sets %>% 
  dplyr::select(c("ID", "Name", "Regions", "trend") 
                | contains("pH6.8_AA") | contains("pH6.4_AA") | contains("pH6.4_CA") | contains("Hypo")) %>%
  dplyr::group_split(Regions, .keep = T)
names(U87.sets) <- names(regions)

# Save the table
sapply(names(U87.sets), function(x){
  openxlsx::write.xlsx(U87.sets[[x]], 
                       file.path(wd, results_dir, date, "tables", 
                                 paste0("U87_GO+pathway_venn_", x, "_", date,".xlsx")))
})

rm(list = c(ls(pattern = "matrix"), "regions", "U87.NES", "U87.FDR"))

# ---------------------------------------------------------------------------- #
palette = c("TGF_BETA_SIGNALING_PATHWAY" = "blue1",
            "EPYTHELIAL_TO_MESENCHYMAL_TRANSITION" = "green3",
            "HYPOXIA" = "brown2")

GSEA.object = list(
  U3017 = HGCC.GSEA$U3017$gsea,
  U3047 = HGCC.GSEA$U3047$gsea,
  U3054 = HGCC.GSEA$U3054$gsea,
  CCLD = CCLD.GSEA$gsea,
  U87_CA = U87.GSEA$`sel_pH647-control_sel`$gsea,
  U87_AA = U87.GSEA$`acu_pH64-control_acu`$gsea
)
GO.object = list(
  U3017 = HGCC.GO$U3017$gsea,
  U3047 = HGCC.GO$U3047$gsea,
  U3054 = HGCC.GO$U3054$gsea,
  CCLD = CCLD.GO$gsea,
  U87_CA = U87.GO$`sel_pH647-control_sel`$gsea,
  U87_AA = U87.GO$`acu_pH64-control_acu`$gsea

)

# ---------------------------------  TGF-B  ---------------------------------- #
tgfb.ranks = list(
  U3017 = which(HGCC.GSEA$U3017$df$ID == "hsa04350"),
  U3047 = which(HGCC.GSEA$U3047$df$ID == "hsa04350"),
  U3054 = which(HGCC.GSEA$U3054$df$ID == "hsa04350"),
  CCLD = which(CCLD.GSEA$df$ID == "hsa04350"),
  U87_CA = which(U87.GSEA$`sel_pH647-control_sel`$df$ID == "hsa04350"),
  U87_AA = which(U87.GSEA$`acu_pH64-control_acu`$df$ID == "hsa04350"))

TGFB.enrichplot <- list()
TGFB.enrichplot$TOTAL_scores  <- do.call(rbind, c(
  lapply(names(tgfb.ranks)[1:4], function(i){
    gsInfo(tgfb.ranks[[i]], object = GSEA.object[[i]]) %>% 
      dplyr::mutate(source = i)}
    )))
TGFB.enrichplot$TOTAL_scores <- TGFB.enrichplot$TOTAL_scores %>% 
  dplyr::mutate(
    Description = gsub(x = TGFB.enrichplot$TOTAL_scores$Description,
                       pattern = "KEGG_", replacement =  ""),
    source = factor(source, levels = names(GSEA.object)[col_order]))

TGFB.enrichplot$TOTAL_table <- getEnrichmentTable(
  .df = rbind(HGCC.GSEA$U3017$df[tgfb.ranks$U3017,] %>% 
                dplyr::mutate(source = "U3017"),
              HGCC.GSEA$U3047$df[tgfb.ranks$U3047,] %>% 
                dplyr::mutate(source = "U3047"),
              HGCC.GSEA$U3054$df[tgfb.ranks$U3054,] %>% 
                dplyr::mutate(source = "U3054"),
              CCLD.GSEA$df[tgfb.ranks$CCLD,] %>% 
                dplyr::mutate(source = "CCLD")),
  .order = col_order, .name = "source")

TGFB.enrichplot$TOTAL_plot <- cowplot::plot_grid(
  plotRunningScore(.df = TGFB.enrichplot$TOTAL_scores, 
                         .x = "x", .y = "runningScore", 
                         .color = "source", .palette = c("U3017" = "pink",
                                                         "U3047" = "salmon",
                                                         "U3054" = "darkred",
                                                         "CCLD" = "steelblue")),
  plotGeneRank(.df = TGFB.enrichplot$TOTAL_scores, 
                     .x = "x", .facet = "source~.",
                     .color = "source", .palette = c("U3017" = "pink",
                                                     "U3047" = "salmon",
                                                     "U3054" = "darkred",
                                                     "CCLD" = "steelblue")),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r", #rel_widths =  c(1,2),
  margins = c(0.5, 0.5, 0.5, 0.5))

results_dir <-  "Results/HGCC"
ggsave(file.path(results_dir, date, "plots", "TOTAL_TGFB_GSEA_enrichment_plot.svg"), bg = "white", device = "svg",
       plot = TGFB.enrichplot$TOTAL_plot, width = 6, height = 4.45, units = "in")

TGFB.enrichplot$U87_scores <- do.call(rbind, c(
  lapply(names(tgfb.ranks)[5:6], function(i){
    gsInfo(tgfb.ranks[[i]], object = GSEA.object[[i]]) %>% 
      dplyr::mutate(source = i)}
  )))
TGFB.enrichplot$U87_scores <- TGFB.enrichplot$U87_scores %>% 
  dplyr::mutate(
    Description = gsub(x = TGFB.enrichplot$U87_scores$Description,
                       pattern = "KEGG_", replacement =  ""),
    source = factor(source, levels = names(GSEA.object)[6:5]))

TGFB.enrichplot$U87_table <- getEnrichmentTable(
  .df = rbind(U87.GSEA$`sel_pH647-control_sel`$df[tgfb.ranks$U87_CA,] %>% 
                dplyr::mutate(source = "U87_CA"),
              U87.GSEA$`acu_pH64-control_acu`$df[tgfb.ranks$U87_AA,] %>% 
                dplyr::mutate(source = "U87_AA")),
  .order = c(2,1), .name = "source")

TGFB.enrichplot$U87_plot <- cowplot::plot_grid(
  plotRunningScore(.df = TGFB.enrichplot$U87_scores, 
                         .x = "x", .y = "runningScore", 
                         .color = "source", .palette = c("U87_CA" = "#F00000",
                                                         "U87_AA" = "salmon")),
  plotGeneRank(.df = TGFB.enrichplot$U87_scores, 
                     .x = "x", .facet = "source~.",
                     .color = "source", .palette = c("U87_CA" = "#F00000",
                                                     "U87_AA" = "salmon")),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r", #rel_widths =  c(1,2),
  margins = c(0.5, 0.5, 0.5, 0.5))

results_dir <-  "Results/U87"
ggsave(file.path(results_dir, date, "plots", "U87_TGFB_GSEA_enrichment_plot.svg"), bg = "white", device = "svg",
       plot = TGFB.enrichplot$U87_plot, width = 6, height = 4.45, units = "in")

# ---------------------------------   EMT   ---------------------------------- #
emt.ranks = list(
  U3017 = which(HGCC.GO$U3017$df$ID == "GO:0001837"),
  U3047 = which(HGCC.GO$U3047$df$ID == "GO:0001837"),
  U3054 = which(HGCC.GO$U3054$df$ID == "GO:0001837"),
  CCLD = which(CCLD.GO$df$ID == "GO:0001837"),
  U87_CA = which(U87.GO$`sel_pH647-control_sel`$df$ID == "GO:0001837"),
  U87_AA = which(U87.GO$`acu_pH64-control_acu`$df$ID == "GO:0001837"))

EMT.enrichplot <- list()
EMT.enrichplot$TOTAL_scores  <- do.call(rbind, c(
  lapply(names(emt.ranks)[1:4], function(i){
    gsInfo(emt.ranks[[i]], object = GO.object[[i]]) %>% 
      dplyr::mutate(source = i)}
  )))
EMT.enrichplot$TOTAL_scores <- EMT.enrichplot$TOTAL_scores %>%
  dplyr::mutate(
    Description = gsub(x = EMT.enrichplot$TOTAL_scores$Description,
                       pattern = "GOBP_", replacement =  ""),
    source = factor(source, levels = names(GO.object)[c(3,1,4,2)]))

EMT.enrichplot$TOTAL_table <- getEnrichmentTable(
  .df = rbind(HGCC.GO$U3017$df[emt.ranks$U3017,] %>% 
                dplyr::mutate(source = "U3017"),
              HGCC.GO$U3047$df[emt.ranks$U3047,] %>% 
                dplyr::mutate(source = "U3047"),
              HGCC.GO$U3054$df[emt.ranks$U3054,] %>% 
                dplyr::mutate(source = "U3054"),
              CCLD.GO$df[emt.ranks$CCLD,] %>% 
                dplyr::mutate(source = "CCLD")),
  .order = col_order, .name = "source")

EMT.enrichplot$TOTAL_plot <- cowplot::plot_grid(
  plotRunningScore(.df = EMT.enrichplot$TOTAL_scores, 
                         .x = "x", .y = "runningScore", 
                         .color = "source", .palette = c("U3017" = "pink",
                                                         "U3047" = "salmon",
                                                         "U3054" = "darkred",
                                                         "CCLD" = "steelblue")),
  plotGeneRank(.df = EMT.enrichplot$TOTAL_scores, 
                     .x = "x", .facet = "source~.",
                     .color = "source", .palette = c("U3017" = "pink",
                                                     "U3047" = "salmon",
                                                     "U3054" = "darkred",
                                                     "CCLD" = "steelblue")),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r", #rel_widths =  c(1,2),
  margins = c(0.5, 0.5, 0.5, 0.5))

results_dir <-  "Results/HGCC"
ggsave(file.path(results_dir, date, "plots", "TOTAL_EMT_GO_enrichment_plot.svg"), bg = "white", device = "svg",
       plot = EMT.enrichplot$TOTAL_plot, width = 6, height = 4.45, units = "in")

EMT.enrichplot$U87_scores <- do.call(rbind, c(
  lapply(names(emt.ranks)[5:6], function(i){
    gsInfo(emt.ranks[[i]], object = GO.object[[i]]) %>% 
      dplyr::mutate(source = i)}
  )))
EMT.enrichplot$U87_scores <- EMT.enrichplot$U87_scores %>%
  dplyr::mutate(
    Description = gsub(x = EMT.enrichplot$U87_scores$Description,
                       pattern = "GOBP_", replacement =  ""),
    source = factor(source, levels = names(GO.object)[6:5]))

EMT.enrichplot$U87_table <- getEnrichmentTable(
  .df = rbind(U87.GO$`sel_pH647-control_sel`$df[emt.ranks$U87_CA,] %>% 
                dplyr::mutate(source = "U87_CA"),
              U87.GO$`acu_pH64-control_acu`$df[emt.ranks$U87_AA,] %>% 
                dplyr::mutate(source = "U87_AA")),
  .order = c(2,1), .name = "source")

EMT.enrichplot$U87_plot <- cowplot::plot_grid(
  plotRunningScore(.df = EMT.enrichplot$U87_scores, 
                         .x = "x", .y = "runningScore", 
                         .color = "source", .palette = c("U87_CA" = "#F00000",
                                                         "U87_AA" = "salmon")),
  plotGeneRank(.df = EMT.enrichplot$U87_scores, 
                     .x = "x", .facet = "source~.",
                     .color = "source", .palette = c("U87_CA" = "#F00000",
                                                     "U87_AA" = "salmon")),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r", #rel_widths =  c(1,2),
  margins = c(0.5, 0.5, 0.5, 0.5))

results_dir <-  "Results/U87"
ggsave(file.path(results_dir, date, "plots", "U87_EMT_GO_enrichment_plot.svg"), bg = "white", device = "svg",
       plot = EMT.enrichplot$U87_plot, width = 6, height = 4.45, units = "in")

# --------------------------------- Hypoxia ---------------------------------- #
hypoxia.ranks = list(
  U3017 = which(HGCC.GO$U3017$df$ID == "M5891"),
  U3047 = which(HGCC.GO$U3047$df$ID == "M5891"),
  U3054 = which(HGCC.GO$U3054$df$ID == "M5891"),
  CCLD = which(CCLD.GO$df$ID == "M5891"),
  U87_CA = which(U87.GO$`sel_pH647-control_sel`$df$ID == "M5891"),
  U87_AA = which(U87.GO$`acu_pH64-control_acu`$df$ID == "M5891"))

Hypoxia.enrichplot <- list()
Hypoxia.enrichplot$TOTAL_scores  <- do.call(rbind, c(
  lapply(names(hypoxia.ranks)[1:4], function(i){
    gsInfo(hypoxia.ranks[[i]], object = GO.object[[i]]) %>% 
      dplyr::mutate(source = i)}
  )))
Hypoxia.enrichplot$TOTAL_scores <- Hypoxia.enrichplot$TOTAL_scores %>%
  dplyr::mutate(
    Description = gsub(x = Hypoxia.enrichplot$TOTAL_scores$Description,
                       pattern = "GOBP_", replacement =  ""),
    source = factor(source, levels = names(GO.object)[col_order]))

Hypoxia.enrichplot$TOTAL_table <- getEnrichmentTable(
  .df = rbind(HGCC.GO$U3017$df[hypoxia.ranks$U3017,] %>% 
                dplyr::mutate(source = "U3017"),
              HGCC.GO$U3047$df[hypoxia.ranks$U3047,] %>% 
                dplyr::mutate(source = "U3047"),
              HGCC.GO$U3054$df[hypoxia.ranks$U3054,] %>% 
                dplyr::mutate(source = "U3054"),
              CCLD.GO$df[hypoxia.ranks$CCLD,] %>% 
                dplyr::mutate(source = "CCLD")),
  .order = col_order, .name = "source")

Hypoxia.enrichplot$TOTAL_plot <- cowplot::plot_grid(
  plotRunningScore(.df = Hypoxia.enrichplot$TOTAL_scores, 
                         .x = "x", .y = "runningScore", 
                         .color = "source", .palette = c("U3017" = "pink",
                                                         "U3047" = "salmon",
                                                         "U3054" = "darkred",
                                                         "CCLD" = "steelblue")),
  plotGeneRank(.df = Hypoxia.enrichplot$TOTAL_scores, 
                     .x = "x", .facet = "source~.",
                     .color = "source", .palette = c("U3017" = "pink",
                                                     "U3047" = "salmon",
                                                     "U3054" = "darkred",
                                                     "CCLD" = "steelblue")),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r", #rel_widths =  c(1,2),
  margins = c(0.5, 0.5, 0.5, 0.5))

results_dir <-  "Results/HGCC"
ggsave(file.path(results_dir, date, "plots", "TOTAL_Hypoxia_GO_enrichment_plot.svg"), bg = "white", device = "svg",
       plot = Hypoxia.enrichplot$TOTAL_plot, width = 6, height = 4.45, units = "in")

Hypoxia.enrichplot$U87_scores <- do.call(rbind, c(
  lapply(names(hypoxia.ranks)[5:6], function(i){
    gsInfo(hypoxia.ranks[[i]], object = GO.object[[i]]) %>% 
      dplyr::mutate(source = i)}
  )))
Hypoxia.enrichplot$U87_scores <- Hypoxia.enrichplot$U87_scores %>%
  dplyr::mutate(
    Description = gsub(x = Hypoxia.enrichplot$U87_scores$Description,
                       pattern = "GOBP_", replacement =  ""),
    source = factor(source, levels = names(GO.object)[5:6]))

Hypoxia.enrichplot$U87_table <- getEnrichmentTable(
  .df = rbind(U87.GO$`sel_pH647-control_sel`$df[hypoxia.ranks$U87_CA,] %>% 
                dplyr::mutate(source = "U87_CA"),
              U87.GO$`acu_pH64-control_acu`$df[hypoxia.ranks$U87_AA,] %>% 
                dplyr::mutate(source = "U87_AA")),
  .order = c(1,2), .name = "source")

Hypoxia.enrichplot$U87_plot <- cowplot::plot_grid(
  plotRunningScore(.df = Hypoxia.enrichplot$U87_scores, 
                         .x = "x", .y = "runningScore", 
                         .color = "source", .palette = c("U87_CA" = "#F00000",
                                                         "U87_AA" = "salmon")),
  plotGeneRank(.df = Hypoxia.enrichplot$U87_scores, 
                     .x = "x", .facet = "source~.",
                     .color = "source", .palette = c("U87_CA" = "#F00000",
                                                     "U87_AA" = "salmon")),
  byrow = T, nrow = 2, ncol = 1, scale = .95, 
  rel_heights = c(1.2,1), axis = "r", #rel_widths =  c(1,2),
  margins = c(0.5, 0.5, 0.5, 0.5))

results_dir <-  "Results/U87"
ggsave(file.path(results_dir, date, "plots", "U87_Hypoxia_GO_enrichment_plot.svg"), bg = "white", device = "svg",
       plot = Hypoxia.enrichplot$U87_plot, width = 6, height = 4.45, units = "in")

################################################################################
# 4. Vulcano visualization of the shared terms and pathways                    #
################################################################################
interest_pathways <- read.csv("data/pathways-of-interest.txt", header = T,
                              sep = "\t", stringsAsFactors = T)

TOTAL.GSEA.vulcano <- list(
  U3017 = rbind(HGCC.GSEA$U3017$sig_df[,c("ID","Name", "NES", "p.adjust")],
                HGCC.GO$U3017$sig_df[,c("ID","Name", "NES","p.adjust")]),
  U3047 = rbind(HGCC.GSEA$U3047$sig_df[,c("ID","Name", "NES","p.adjust")],
                HGCC.GO$U3047$sig_df[,c("ID","Name", "NES","p.adjust")]),
  U3054 = rbind(HGCC.GSEA$U3054$sig_df[,c("ID","Name", "NES","p.adjust")],
                HGCC.GO$U3054$sig_df[,c("ID","Name", "NES","p.adjust")]),
  CCLD = rbind(CCLD.GSEA$sig_df[,c("ID","Name", "NES","p.adjust")],
               CCLD.GO$sig_df[,c("ID","Name", "NES","p.adjust")]))

TOTAL.GSEA.vulcano <- lapply(TOTAL.GSEA.vulcano, function(x){
  return(x %>% 
           dplyr::mutate(p.adjust = ifelse(p.adjust < 1e-35, yes = 1e-35, no = p.adjust))
           )
})

### Extracellular matrix organization
(TOTAL.ECM.vulcano <- lapply(TOTAL.GSEA.vulcano, function(x){
  plotClusters(.df = x, 
               .pathways = interest_pathways %>% filter(Category == "ECM"))
}))

results_dir <- "Results/HGCC"
sapply(names(TOTAL.ECM.vulcano), function(x){
  ggsave(file.path(results_dir,date,"plots",
                   paste(x,"ECM_GSEA_Vulcano_plot.svg",sep = "_")),
         plot = TOTAL.ECM.vulcano[[x]], bg = "white", device = "svg",
         width = 12, height = 8,units = "in")
})

TOTAL.ECM.vulcano$HGCC <- cowplot::plot_grid(
  TOTAL.ECM.vulcano$U3017, TOTAL.ECM.vulcano$U3047, TOTAL.ECM.vulcano$U3054,
  byrow = T, nrow = 1, ncol = 3, scale = .9, 
  margins = c(0.5, 0.5, 0.5, 0.5))

ggsave(file.path(results_dir,date,"plots",
                 paste("HGCC_combined_ECM_GSEA_Vulcano_plot.svg",sep = "_")),
       plot = TOTAL.ECM.vulcano$HGCC, bg = "white", device = "svg",
       width = 20, height = 8, units = "in")

### Hypoxia
(TOTAL.HYPO.vulcano <- lapply(TOTAL.GSEA.vulcano, function(x){
  plotClusters(.df = x, 
               .pathways = interest_pathways %>% filter(Category == "HYPOXIA"))
}))

sapply(names(TOTAL.HYPO.vulcano), function(x){
  ggsave(file.path(results_dir,date,"plots",
                   paste(x,"Hypoxia_GSEA_Vulcano_plot.svg",sep = "_")),
         plot = TOTAL.HYPO.vulcano[[x]], bg = "white", device = "svg",
         width = 12, height = 8,units = "in")
})

TOTAL.HYPO.vulcano$HGCC <- cowplot::plot_grid(
  TOTAL.HYPO.vulcano$U3017, TOTAL.HYPO.vulcano$U3047, TOTAL.HYPO.vulcano$U3054,
  byrow = T, nrow = 1, ncol = 3, scale = .9, 
  margins = c(0.5, 0.5, 0.5, 0.5))

ggsave(file.path(results_dir,date,"plots",
                 paste("HGCC_combined_Hypoxia_GSEA_Vulcano_plot.svg",sep = "_")),
       plot = TOTAL.HYPO.vulcano$HGCC, bg = "white", device = "svg",
       width = 20, height = 8, units = "in")

### TGF-beta signaling pathway
(TOTAL.TGFb.vulcano <- lapply(TOTAL.GSEA.vulcano, function(x){
  plotClusters(.df = x, 
               .pathways = interest_pathways %>% filter(Category == "TGFb"))
}))

sapply(names(TOTAL.TGFb.vulcano), function(x){
  ggsave(file.path(results_dir,date,"plots",
                   paste(x,"TGFb_GSEA_Vulcano_plot.svg",sep = "_")),
         plot = TOTAL.TGFb.vulcano[[x]], bg = "white", device = "svg",
         width = 12, height = 8,units = "in")
})

TOTAL.TGFb.vulcano$HGCC <- cowplot::plot_grid(
  TOTAL.TGFb.vulcano$U3017, TOTAL.TGFb.vulcano$U3047, TOTAL.TGFb.vulcano$U3054,
  byrow = T, nrow = 1, ncol = 3, scale = .9, 
  margins = c(0.5, 0.5, 0.5, 0.5))

ggsave(file.path(results_dir,date,"plots",
                 paste("HGCC_combined_TGFb_GSEA_Vulcano_plot.svg",sep = "_")),
       plot = TOTAL.TGFb.vulcano$HGCC, bg = "white", device = "svg",
       width = 20, height = 8, units = "in")

# ---------------------------------   U87   ---------------------------------- #
U87.GSEA.vulcano <- list(
  U87_CA64 = rbind(U87.GSEA$`sel_pH647-control_sel`$sig_df[,c("ID","Name", "NES","p.adjust")],
                 U87.GO$`sel_pH647-control_sel`$sig_df[,c("ID","Name", "NES","p.adjust")]),
  U87_AA64 = rbind(U87.GSEA$`acu_pH64-control_acu`$sig_df[,c("ID","Name", "NES","p.adjust")],
                 U87.GO$`acu_pH64-control_acu`$sig_df[,c("ID","Name", "NES","p.adjust")]),
  U87_AA68 = rbind(U87.GSEA$`acu_pH68-control_acu`$sig_df[,c("ID","Name", "NES","p.adjust")],
                 U87.GO$`acu_pH68-control_acu`$sig_df[,c("ID","Name", "NES","p.adjust")]),
  U87_Hypoxia = rbind(U87.GSEA$`hypoxia-control_nox`$sig_df[,c("ID","Name", "NES","p.adjust")],
                 U87.GO$`hypoxia-control_nox`$sig_df[,c("ID","Name", "NES","p.adjust")])
)

U87.GSEA.vulcano <- lapply(U87.GSEA.vulcano, function(x){
  return(x %>% 
           dplyr::mutate(p.adjust = ifelse(p.adjust < 1e-35, yes = 1e-35, no = p.adjust))
  )
})

### Extracellular matrix organization
(U87.ECM.vulcano <- lapply(U87.GSEA.vulcano, function(x){
  plotClusters(.df = x, 
               .pathways = interest_pathways %>% filter(Category == "ECM"))
}))

results_dir <- "Results/U87"
sapply(names(U87.ECM.vulcano), function(x){
  ggsave(file.path(results_dir,date,"plots",
                   paste(x,"ECM_GSEA_Vulcano_plot.svg",sep = "_")),
         plot = U87.ECM.vulcano[[x]], bg = "white", device = "svg",
         width = 12, height = 8,units = "in")
})

### Hypoxia
(U87.HYPO.vulcano <- lapply(U87.GSEA.vulcano, function(x){
  plotClusters(.df = x, 
               .pathways = interest_pathways %>% filter(Category == "HYPOXIA"))
}))

sapply(names(U87.HYPO.vulcano), function(x){
  ggsave(file.path(results_dir,date,"plots",
                   paste(x,"Hypoxia_GSEA_Vulcano_plot.svg",sep = "_")),
         plot = U87.HYPO.vulcano[[x]], bg = "white", device = "svg",
         width = 12, height = 8,units = "in")
})

### TGF-beta signaling pathway
(U87.TGFb.vulcano <- lapply(U87.GSEA.vulcano, function(x){
  plotClusters(.df = x, 
               .pathways = interest_pathways %>% filter(Category == "TGFb"))
}))

sapply(names(U87.TGFb.vulcano), function(x){
  ggsave(file.path(results_dir,date,"plots",
                   paste(x,"TGFb_GSEA_Vulcano_plot.svg",sep = "_")),
         plot = U87.TGFb.vulcano[[x]], bg = "white", device = "svg",
         width = 12, height = 8,units = "in")
})
