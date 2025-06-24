################################################################################
# Created: 2024 10 18 ; Last Modified: 2025 06 24 ; MH                         #
################################################################################
# 1.) Comparing the Uppsala cells and Lipid loaded cells
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

interest_genes <- read.csv("../data/gene-signature-new.txt", header = T,
                           sep = "\t", stringsAsFactors = T)
# interest_genes:
# "CA9","VEGFA","CA12","BGN","CSPG4","VCAN","NCAN","CHPF","CSGALNACT1","CHSY1","UST",
# "SULF2","COL6A1","FN1","THBS4","GPC4","FASN","PPARD","PPARGC1A","VLDLR","HILPDA"
TOTAL.heatmap <- list()
TOTAL.heatmap$data <- interest_genes %>% 
  dplyr::inner_join(., CCLD.df[,1:2], by = c("SYMBOL" = "Symbol")) %>%
  dplyr::inner_join(., HGCC.deg$U3017[,c(2,4)], by = c("SYMBOL" = "Symbol")) %>% 
  dplyr::inner_join(., HGCC.deg$U3047[,c(2,4)], by = c("SYMBOL" = "Symbol")) %>%
  dplyr::inner_join(., HGCC.deg$U3054[,c(2,4)], by = c("SYMBOL" = "Symbol")) %>% 
  dplyr::rename_at(vars(contains("log2")), ~c("LDvsnoLD", "U3017", "U3047", "U3054")) %>%
  dplyr::mutate(Category = factor(Category, 
                                      levels = c("Acidosis/Hypoxia", "CSPG core", "CSPG Biosynthesis",
                                                 "ECM remodelling", "Lipid metabolism"),
                                      labels = c("Acidosis/Hypoxia", "CSPG core", "CSPG Biosynthesis",
                                                 "ECM remodelling", "Lipid metabolism"))) %>% 
  dplyr::select(SYMBOL, Category, U3017, U3047, U3054, LDvsnoLD)

TOTAL.heatmap$matrix <- as.matrix(TOTAL.heatmap$data[,c(3:6)])

row_order <- match(c("CA9","VEGFA","CA12","BGN","CSPG4","VCAN","NCAN",
                     "CHPF","CSGALNACT1","CHSY1","UST","SULF2","COL6A1","FN1",
                     "THBS4","GPC4","FASN","PPARD","PPARGC1A","HILPDA","VLDLR"),
                   TOTAL.heatmap$data$SYMBOL)
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


row.names(TOTAL.heatmap$matrix) <- TOTAL.heatmap$data$SYMBOL 
TOTAL.heatmap$matrix <- TOTAL.heatmap$matrix[row_order,col_order]#row_order

cat_annot = rowAnnotation(
  Category = TOTAL.heatmap$data$Category[row_order],
  col = list(Category = c("Acidosis/Hypoxia" = "darkgreen",
                          "CSPG core" = "turquoise",
                          "CSPG Biosynthesis" = "#FF08FF",
                          "ECM remodelling" = "#8D00FA",
                          "Lipid metabolism" = "white")
            
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
                              right_annotation = cat_annot, 
                              heatmap_legend_param = list(
                                title = "Log2 (Fold change)",
                                direction = "vertical"),
                              border_gp = gpar(col = "black", lty = 2)))

svg(file.path("Results","signature", "category_heatmap_new_signature.svg"), width = 12, height = 8)
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

