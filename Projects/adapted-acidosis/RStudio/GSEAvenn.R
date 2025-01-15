################################################################################
# 1.) The comparison of the 3 Uppsala cell lines HGCC)                         #
################################################################################
# 1. Load the GSEA results for the 3 cell lines
tmp <- list(HGCC.GSEA$U3017$sig_df[,c("ID", "NES","p.adjust")],
            HGCC.GSEA$U3047$sig_df[,c("ID", "NES", "p.adjust")],
            HGCC.GSEA$U3054$sig_df[,c("ID", "NES", "p.adjust")])
names(tmp) <- names(HGCC.GSEA)

# Extract the number of genes in each region of the Venn diagram
sets <- GOVenn(HGCC.GSEA$U3017$sig_df[,c("ID", "NES","p.adjust")],
               HGCC.GSEA$U3047$sig_df[,c("ID", "NES","p.adjust")],
               HGCC.GSEA$U3054$sig_df[,c("ID", "NES","p.adjust")],
               plot = F)$table
sets <- lapply(sets, rownames_to_column, var = "ID")
sets[["A_only"]] <- sets[["A_only"]] %>% dplyr::rename(logFC_A = logFC)
sets[["B_only"]] <- sets[["B_only"]] %>% dplyr::rename(logFC_B = logFC)
sets[["C_only"]] <- sets[["C_only"]] %>% dplyr::rename(logFC_C = logFC)

# Create a data frame for the Venn diagram
table <- do.call(rbind.fill,rev(sets))
table <- table %>% dplyr::distinct(ID, .keep_all = T)

# Create an upset plot base: change log2FC values to logical (TRUE/FALSE)
table <- table %>%
  dplyr::mutate("U3017" = ifelse(is.na(logFC_A), FALSE, TRUE)) %>%
  dplyr::mutate("U3047" = ifelse(is.na(logFC_B), FALSE, TRUE)) %>%
  dplyr::mutate("U3054" = ifelse(is.na(logFC_C), FALSE, TRUE)) %>%
  dplyr::relocate(where(is.logical), .before = where(is.character))
colnames(table)[5:7] = c("NES_U3017", "NES_U3047", "NES_U3054")

# get the genes in each of the regions
arranged <- arrange_venn(table, c("U3017","U3047","U3054"), radius = 2, max_iterations = 500,
                         extract_regions = T)
# arranged_df <- arrange_venn(table, c("U3017","U3047","U3054"),
#                             extract_regions = F, radius = 2, max_iterations = 1000)
set.seed(123)
# randomize the position of the genes in the regions
xy = rbind(
  calculateCircle(x = 0.00, y = 0.2886, r = 0.1, # 3-way intersection
                  noiseFun = function(x) (x + rnorm(1,0,0.1)), # add noise
                  steps = 54,randomDist = T, randomFun = rnorm), # 579 genes
  calculateEllipse(x = -0.65, y = -0.0866, a = 0.2, b = .2, angle = -240, 
                   noiseFun = function(x) (x + rnorm(1,0,0.1)),
                   steps = 22, randomDist = T, randomFun = rnorm),
  calculateEllipse(x = 0.65, y = -0.0866, a = 0.2, b = .2, angle = -120, 
                   noiseFun = function(x) (x + rnorm(1,0,0.1)),
                   steps = 68, randomDist = T, randomFun = rnorm),
  calculateEllipse(x = 0, y = 1.0392, a = 0.1, b = .2, angle = 0, 
                   noiseFun = function(x) (x + rnorm(1,0,0.1)),
                   steps = 17, randomDist = T, randomFun = rnorm),
  calculateCircle(x = 0.00, y = -1.2124, r = .3, 
                  noiseFun = function(x) (x + rnorm(1,0,0.2)),
                  steps = 131, randomDist = T, randomFun = rnorm),
  calculateCircle(x = 1.30, y = 1.0392, r = .3,  
                  noiseFun = function(x) (x + rnorm(1,0,0.2)),
                  steps = 29,randomDist = T, randomFun = rnorm),
  calculateCircle(x = -1.30, y = 1.0392, r = .3, 
                  noiseFun = function(x) (x + rnorm(1,0,0.2)),
                  steps = 62,randomDist = T, randomFun = rnorm)
)

table <- table %>%
  dplyr::mutate(region = dplyr::case_when(
    `U3017` & `U3047` & `U3054` ~ "U3017-U3047-U3054",
    `U3017` & !`U3047` & `U3054` ~ "U3017-U3054",
    !`U3017` & `U3047` & `U3054` ~ "U3047-U3054",
    `U3017` & `U3047` & !`U3054` ~ "U3017-U3047",
    !`U3017` & !`U3047` & `U3054` ~ "U3054",
    !`U3017` & `U3047` & !`U3054` ~ "U3047",
    `U3017` & !`U3047` & !`U3054` ~ "U3017"
  ), # relabel the regions
  region = as.factor(region),
  x = xy[,1],
  y = xy[,2],
  )

# 2. Create a Venn diagram of the GSEA results
HGCC_GSEA_venn <- list()
(HGCC_GSEA_venn$plot <- ggplot(table)
 # --- BASE PLOT & LEGEND--- #
 # create the base plot with 1:1 aspect ratio and clear whit background
 + coord_fixed()
 + theme_void()
 # Add invisible gene sets to be able to add the legend to the plot
 + geom_point(aes(x = x, y = y, colour = Trend), 
              alpha = 0.01, size = 1)
 + scale_color_manual(
   labels = c("Opposing regulation","Inhibited","Activated"),
   values = c("Change" = "orange", "DOWN" = "blue", "UP" = "darkred"),
   guide = guide_legend(override.aes = c(alpha = 1, size = 5)),
   aesthetics = "colour",
 ) 
 # --- SETS: BORDERS, FILL, LABELS --- #
 # Add the set borders of the Venn diagram
 + geom_venn_region(table, sets = c("U3017","U3047","U3054"), 
                    alpha = 0.3, show.legend = F) 
 + geom_venn_circle(table, sets = c("U3017","U3047","U3054"),
                    size = .5)
 # Fill every set white
 + scale_fill_venn_mix(table, sets =  c("U3017","U3047","U3054"),
                       guide='none', highlight = c("U3017-U3047-U3054"),
                       inactive_color='NA', active_color = 'NA')
 # Add the names of the sets
 + geom_venn_label_set(table, outwards_adjust = 2,
                       sets = c("U3017","U3047","U3054"),
                       fill = alpha("black", .35), 
                       aes(label =  c("U3054","U3047","U3017")))
 # --- DECORATE SETS; GENE EXPRESSION TRENDS --- #
 # add labels with numbers of incoherently activated genes
 + geom_venn_label_region(
   table %>% dplyr::filter(Trend == "Change"),
   sets = c("U3017","U3047","U3054"),
   aes(label=size ), fill = alpha("orange", .35),
   outwards_adjust=1.25, position=position_nudge(y=0.15))
 # add the numbers of up-regulated genes
 + geom_venn_label_region(
   table %>% dplyr::filter(Trend == "UP"),
   sets = c("U3017","U3047","U3054"),
   aes(label=size ), fill = alpha("darkred", .35),
   outwards_adjust=1.25, position=position_nudge(y=-0.15, x = 0.2))
 # add the numbers of down-regulated genes
 + geom_venn_label_region(
   table %>% dplyr::filter(Trend == "DOWN"),
   sets = c("U3017","U3047","U3054"),
   aes(label=size), fill = alpha("blue", .35),
   outwards_adjust=1.25, position=position_nudge(y=-0.15, x = -0.2)) 
 # --- SET THE THEME --- #
 + theme(
   legend.title = element_text(size = 12),
   legend.text = element_text(size = 12)))

# Save the plot
results_dir <- "Results/HGCC"
ggsave(file.path(wd, results_dir, "HGCC_GSEA_venn.png"), bg = "white", 
       plot = HGCC_GSEA_venn$plot, width = 12, height = 8)

# 3. Extract the genes in each region of the Venn diagram
tmp <- merge.rec(tmp, by = "ID", all = T, suffixes = c("", ""))
colnames(tmp) <- c("ID", "NES (U3017)", "p.adjust (U3017)", 
                   "NES (U3047)", "p.adjust (U3047)", "NES (U3054)", "p.adjust (U3054)")

HGCC_GSEA_venn$table <- table %>% 
  dplyr::select(ID, Trend, region) %>% 
  dplyr::inner_join(tmp, by = "ID") %>% 
  dplyr::mutate(Trend = dplyr::case_when(
    Trend == "Change" ~ "Opposing regulation",
    Trend == "UP" ~ "Activated",
    Trend == "DOWN" ~ "Inhibited"
  ),
  Description = pathways$gs_description[match(ID, pathways$gs_exact_source)]) %>% 
  dplyr::relocate(Description, .after = ID)

HGCC_GSEA_venn$table <- HGCC_GSEA_venn$table %>% 
  dplyr::group_split(region, .keep = T)

HGCC_GSEA_venn$table <- setNames(HGCC_GSEA_venn$table, c("U3017","U3017-U3047","U3017-U3047-U3054", "U3017-U3054",
                                 "U3047","U3047-U3054","U3054"))

# Save the table
sapply(names(HGCC_GSEA_venn$table), function(x){
  openxlsx::write.xlsx(HGCC_GSEA_venn$table[[x]], 
                       file.path(results_dir, paste0("HGCC_GSEA_venn_", x, ".xlsx")))
})

################################################################################
# 2.) Comparing the Uppsala cells, U87 acidosis cels, and Lipid loaded cells   #
################################################################################

# GSEA on the GO terms from MSigDB
sets <- list(rbind(HGCC.GSEA$U3017$sig_df[,c("ID","Name", "NES")],
                   HGCC.GO$U3017$sig_df[,c("ID","Name", "NES")]),
             rbind(HGCC.GSEA$U3047$sig_df[,c("ID","Name", "NES")],
                   HGCC.GO$U3047$sig_df[,c("ID","Name", "NES")]),
             rbind(HGCC.GSEA$U3054$sig_df[,c("ID","Name", "NES")],
                   HGCC.GO$U3054$sig_df[,c("ID","Name", "NES")]),
             rbind(CCLD.GSEA$sig_df[,c("ID","Name", "NES")],
                   CCLD.GO$sig_df[,c("ID","Name", "NES")]))
# Merge the data frames
table <- merge.rec(sets, by = c("ID","Name"), suffixes = c("",""), all = T) %>% 
  setNames(., c("ID","Name", c("U3017","U3047","U3054","CCLD")))

# Extract the allocation of each gene
regions <- table %>%
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
table <- dplyr::inner_join(table, regions, by = c("Name", "ID"))

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
point.matrix <- dplyr::inner_join(table, label.matrix, 
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

trend.matrix <- point.matrix %>% 
  dplyr::group_by(Regions, trend) %>% 
  dplyr::summarise(Size = n()) %>% 
  dplyr::inner_join(label.matrix[,-2], by = "Regions") %>% 
  tidyr::pivot_wider(names_from = trend, values_from = Size) %>% 
  dplyr::mutate_all(~ifelse(is.na(.),0,.))

  
# Create the Venn diagram
TOTAL_GSEA_venn <- list()
(TOTAL_GSEA_venn$plot <- ggplot() +
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
       bg = "white", plot = TOTAL_GSEA_venn$plot, width = 12, height = 8)

table <- table %>% dplyr::group_split(Regions, .keep = T)
names(table) <- names(regions)

# Save the table
sapply(names(table), function(x){
  openxlsx::write.xlsx(table[[x]], 
                       file.path(results_dir,"..", 
                                 paste0("TOTAL_GO+pathway_venn_", x,"_", Sys.Date(),".xlsx")))
})

################################################################################
# 3. Visualize selected enrichment scores                                      #
################################################################################

library(enrichplot)
palette = c("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION" = "#8338EC",
            "REACTOME_GLYCOSAMINOGLYCAN_METABOLISM" = "#118AB2",
            "REACTOME_ECM_PROTEOGLYCANS" = "#FF006E",
            "REACTOME_CHONDROITIN_SULFATE_DERMATAN_SULFATE_METABOLISM" = "#FFBE0B")
ranks = which(CCLD.GSEA$df$Name %in% gsub(x = names(palette), "REACTOME_", ""))

gsInfo <- function (object, geneSetID) {
  geneList <- object@geneList
  if (is.numeric(geneSetID)) 
    geneSetID <- object@result[geneSetID, "ID"]
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify = TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  if (length(object@gene2Symbol) == 0) {
    df$gene <- names(geneList)
  }
  else {
    df$gene <- object@gene2Symbol[names(geneList)]
  }
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}
gseaScores <- function (geneList, geneSet, exponent = 1, fortify = FALSE) 
{
  geneSet <- intersect(geneSet, names(geneList))
  N <- length(geneList)
  Nh <- length(geneSet)
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit/NR)
  Pmiss[!hits] <- 1/(N - Nh)
  Pmiss <- cumsum(Pmiss)
  runningES <- Phit - Pmiss
  max.ES <- max(runningES)
  min.ES <- min(runningES)
  if (abs(max.ES) > abs(min.ES)) {
    ES <- max.ES
  }
  else {
    ES <- min.ES
  }
  df <- data.frame(x = seq_along(runningES), runningScore = runningES, 
                   position = as.integer(hits))
  if (fortify == TRUE) {
    return(df)
  }
  df$gene = names(geneList)
  res <- list(ES = ES, runningES = df)
  return(res)
}

runningScores  <- do.call(rbind, lapply(ranks, gsInfo, object = CCLD.GSEA$gsea))
runningScores$Description <- factor(runningScores$Description, levels = names(palette))

p1 <- ggplot(runningScores, aes(x = x)) + 
  xlab(element_blank()) + 
  ylab("Enrichment score (ES)") + theme_bw(base_size = 14) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0.01, 0.01), breaks = c(-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  geom_point(aes(y = runningScore, color = Description),show.legend = F,
             size = 2, data = subset(runningScores, position == 1)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 14)) + 
  scale_color_manual(values = palette)

p2 <- ggplot(runningScores, aes(x = x)) +
  geom_linerange(aes(ymin = ymin, ymax = ymax, color = Description), show.legend = F, 
                 size = 1, data = subset(runningScores, position == 1)) +
  theme_bw(base_size = 14) + 
  xlab("Rank in ordered gene list") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  scale_color_manual(values = palette) +
  facet_grid(Description~., scales = "free_y") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        strip.text = element_blank(), axis.text.x = element_text(colour = "black", size = 14))

mytable <- CCLD.GSEA$sig_df[ranks,] %>% 
  dplyr::arrange(p.adjust) %>% 
  dplyr::mutate(
    NES = round(NES, 4),
    FDR = round(p.adjust, 4)) %>%
  dplyr::select(c(5,15))

table_plot <- tableGrob(mytable, 
                        theme = ttheme_minimal(base_size = 14,
                                               core = list(fg_params = list(hjust = 0.5, x = 0.5,
                                                                col = palette)),
                                               rowhead = list(fg_params = list(hjust = 0, x = 0,
                                                              col = palette[c(4,1,2,3)])))) 

CCLD.enrichplot <- cowplot::plot_grid("", "", p1, "", p2, table_plot, byrow = T,
                   nrow = 3, ncol = 2, 
                   rel_heights = c(1, 2,1), rel_widths =  c(1,1.5))


ggsave(file.path("..", results_dir, "CCLD_GSEA.svg"), bg = "white", device = "svg",
       plot = CCLD.enrichplot, width = 20, height = 8, units = "in")

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
