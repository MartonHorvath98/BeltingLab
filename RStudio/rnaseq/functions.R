# ------------------- Data processing functions --------------------------------
# Recursively merges a list of data frames into a single data frame
merge.rec <- function(.list, ...){
  if(length(.list)==1) return(.list[[1]])
  Recall(c(list(merge(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
}

# Function to bind GO and KEGG results into a single data frame
bind_df <- function(go, kegg){
  go <- lapply(go, function(x){
    x %>% 
      dplyr::select(ONTOLOGY, ID, Description, GeneRatio, pvalue, geneID) %>%
      dplyr::rename("category" = ONTOLOGY)
  })
  
  kegg <- kegg %>% 
    dplyr::select(ID, Description,GeneRatio, pvalue, geneID) %>% 
    dplyr::mutate(category = "KEGG")
  tmp <- do.call(rbind, append(go, list("KEGG" = kegg)))
  return(tmp)
}

# ----------------- Differential expression analysis functions -----------------
# Function to perform differential expression analysis using DESeq2
make_deseq <- function(matrix, coldata, design){
  # Create a formula for the design matrix
  .design <- reformulate(design)
  # Create a DESeqDataSet object
  deseq <- DESeqDataSetFromMatrix(matrix, colData = coldata, design = .design)
  # Perform differential expression analysis
  dds <- DESeq(deseq)
  # Remove features with low counts, <1 read per million in at least 3 samples
  subset <- rowSums(cpm(counts(dds), normalized = T) >1 ) >= 3
  dds <- dds[subset,]
  
  # Perform variance stabilizing transformation
  rld <- rlog(dds, blind = F)
  
  return(list(deseq = deseq, dds = dds, rld = rld))
  rm(list(deseq, dds, subset, rld))
}

# Function to extract a result table from the DESeq analysis
get_results <- function(dds, sig_log2FC, sig_pval, contrast = NULL, name = NULL){
  # Extract results table
  tmp <- results(dds, contrast = contrast, name = name,
                 independentFiltering = T, pAdjustMethod = "BH", alpha = 0.05)
  # Perform multiple testing correction
  fdr <- fdrtool(tmp$stat, statistic = "normal")
  tmp$padj <- p.adjust(fdr$pval, method = "BH")
  # Create a data frame from the results table
  df <- data.frame(tmp, 
                   ensemblID = rownames(assay(dds)), 
                   row.names = rownames(assay(dds))) %>% 
    # Annotate genes with gene symbols and entrez IDs
    dplyr::mutate(
      geneID = mapIds(org.Hs.eg.db, ensemblID, 
                      keytype = "ENSEMBL", column = "SYMBOL",
                      multiVals = "first"),
      entrezID = mapIds(org.Hs.eg.db, ensemblID, 
                        keytype = "ENSEMBL", column = "ENTREZID",
                        multiVals = "first"),
      # Add a column for significance using the specified thresholds
      significance = dplyr::case_when(
        abs(log2FoldChange) > sig_log2FC & padj > sig_pval ~ 'log2FoldChange',
        abs(log2FoldChange) < sig_log2FC & padj < sig_pval ~ 'Log10P',
        log2FoldChange < (-1)*sig_log2FC & padj < sig_pval ~ 'Signif. down-regulated',
        log2FoldChange > sig_log2FC & padj < sig_pval ~ 'Signif. up-regulated',
        T ~ 'NS')) %>%
    dplyr::relocate(c("ensemblID","geneID","entrezID"), .after = everything()) %>% 
    dplyr::filter(complete.cases(.))
  
  # Create a data frame with only significant results
  sig_df <- df %>%
    dplyr::filter(significance == 'Signif. down-regulated' | 
                    significance == 'Signif. up-regulated') %>%
    # Remove unnamed genes
    dplyr::filter(!is.na(geneID))
  
  return(list(results = tmp, df = df, sig_df = sig_df))
  rm(list(tmp, fdr, df, significance, sig_df))
}

# ------------------- Data visualisation functions -----------------------------
# Function to create a PCA plot with biplot
make_pca <- function(rld, group, labs, cols){
  # Create a scaled and mean normalized PCA object
  tmp <- prcomp(t(assay(rld)), center = TRUE, scale. = TRUE)
  
  # Extract the loadings and select the top 5 contributing genes for each PC
  loadings <- tmp$rotation
  top_PC1 <- names(sort(abs(loadings[,1]), decreasing = T)[1:5])
  top_PC2 <- names(sort(abs(loadings[,2]), decreasing = T)[1:5])
  top_genes <- unique(c(top_PC1, top_PC2))
  
  # Create a data frame for plotting the gene arrows
  arrows <- subset(tmp$rotation, 
                   rownames(tmp$rotation) %in% top_genes,
                   c("PC1","PC2"))
  arrows <- as.data.frame(arrows)
  # Calculate the angle of the arrows
  arrows$angle <- with(arrows, (180/pi) * atan(PC1/PC2))
  # Calculate the x and y coordinates for the arrowheads
  arrows$xvar <- with(arrows, 2 * cos(angle))
  arrows$yvar <- with(arrows, 2 * sin(angle))
  # Annotate the arrows with gene symbols
  arrows$geneID <- AnnotationDbi::mapIds(org.Hs.eg.db, rownames(arrows), 
                                         keytype = "ENSEMBL", column = "SYMBOL",
                                         multiVals = "first")
  # Create a ggplot object
  arrow_style <- arrow(length = unit(1/2, "picas"), type = "closed")
  
  # Create the biplot
  plot <- ggbiplot(pcobj = tmp, choices = 1:2, scale = 1, circle = T,
                   groups = rld@colData[[group]], var.axes = F) +
    geom_point(size = 5, aes(color = rld@colData[[group]])) +
    # add confidence ellipses
    stat_ellipse(geom = "polygon", 
                 aes(color = rld@colData[[group]], fill = rld@colData[[group]]), 
                 type = "norm", level = 0.68, 
                 alpha = .3, linetype = 2, show.legend = F) +
    scale_color_manual(name = "Conditions",
                       labels = labs,
                       values = cols) +
    scale_fill_manual(name = "Conditions",
                      labels = labs,
                      values = cols) +
    guides(color=guide_legend(ncol=2),
           fill=guide_legend(ncol=2)) +
    theme(axis.text = element_text(size = 12, colour = "darkgrey"),
          axis.title = element_text(size = 12, face = "bold", colour = "black"),
          legend.text = element_text(size = 12, colour = "darkgrey"),
          legend.title = element_text(size = 12, face = "bold", colour = "black"),
          legend.position = "bottom", legend.direction = "horizontal")
  # Add the gene arrows to the plot
  plot <- plot + 
    geom_segment(data = arrows, 
                 aes(x = 0, y = 0, xend = xvar, yend = yvar), show.legend = F,
                 arrow = arrow_style) + 
    geom_text_repel(aes(label = rld@colData@rownames, color = rld@colData[[group]]), 
                    size = 3, show.legend = F) + 
    geom_label_repel(data = arrows, aes(x = xvar, y = yvar, label = geneID),
                     color = "black", size = 3, arrow = NULL, vjust = 0.5, hjust = 0.5)
  
  return(plot)
  rm(list(tmp, plot))
}

# Function to create vulcan-plot visualisation for the differential expression 
# analysis results
make_vulcanplot <- function(total, sig){
  # color pallete for NS genes, genes that pass p-value or log2FC filters, and
  # significantly down- or up-regulated genes, respectively
  colourPalette <- c('#c1c1c1', '#363636','#222222', '#000f64','#841f27')
  
  # Create a ggplot object
  return((ggplot(data = na.omit(total), 
                 aes(x = log2FoldChange, 
                     y = -log10(padj), 
                     colour = significance)) 
          + geom_point(mapping = aes(), inherit.aes = T, size = 2.5) 
          + scale_color_manual(values = c(
            "NS" = "#c1c1c1",
            "Log10P" = '#363636',
            "Log2FoldChange" = '#767676',
            "Signif. up-regulated" = '#841f27',
            "Signif. down-regulated" = '#000f64'
          ))
          + labs(x = expression(paste(log[2], 'FoldChange')),
                 y = expression(paste(log[10], italic('FDR')))) 
          + theme(axis.title = element_text(size = 14), 
                  axis.text = element_text(size = 14), 
                  legend.position = 'none') 
          + scale_x_continuous(expand = expansion(0.2))
          # Visualize log2FC threshold
          + geom_vline(xintercept = c(-1.5, 1.5), 
                       linetype = 'dotted', size = 1) 
          # Visualize adjusted p-value threshold
          + geom_hline(yintercept = -log10(0.05), 
                       linetype = 'dotted', size = 1)
          # Add labels for significantly up- or down-regulated genes
          + geom_text(data = na.omit(sig), hjust = 0, vjust = 1.5, 
                      colour = 'black', position = 'identity', 
                      show.legend = F, check_overlap = T,
                      label = na.omit(sig)[,"geneID"])))
}

# ------------- Functions for overrepresentation analyses ----------------------
# Function to extract gene lists for KEGG pathway enrichment analysis
get_genelist <- function(list){
  # Extract the background gene list of every expressed gene
  background <- list$df %>%
    dplyr::arrange(desc(log2FoldChange)) %>%
    dplyr::distinct(entrezID, .keep_all = T) %>%
    dplyr::pull("log2FoldChange", name ="entrezID")
  # Extract the gene list of interest of DEGs
  interest <- list$sig_df %>%
    dplyr::pull("log2FoldChange", name ="entrezID")
  interest <- sort(interest, decreasing = T) # sort the gene list
  
  return(list(background = background, interest = interest))
}

# Function to perform KEGG pathway enrichment analysis
kegg_results <- function(x, y){
  # Perform KEGG pathway enrichment analysis
  kegg <- enrichKEGG(names(x), # gene set of interest
                     organism = 'hsa', # human pathway database
                     keyType = 'kegg',
                     universe = names(y), # background gene set
                     pvalueCutoff = 0.05) # return only significant results
  
  # Transform the results into a data frame
  return(setReadable(kegg, org.Hs.eg.db,"ENTREZID") %>%
           as.data.frame(.) %>% 
           dplyr::mutate(
             # Change gene ratio to numeric format
             GeneRatio = sapply(stringr::str_split(GeneRatio, "/"), 
                                function(y) as.numeric(y[1])/as.numeric(y[2])),
             # Change background ratio to numeric format
             BgRatio = sapply(stringr::str_split(BgRatio, "/"), 
                              function(y) as.numeric(y[1])/as.numeric(y[2])),
             Count = as.numeric(Count)))
}

# Function to perform GO enrichment analysis
go_results <- function(x, y, type = c("ALL","BP","MF","CC")){
  # Perform GO enrichment analysis
  go <- enrichGO(names(x), # gene set of interest
                 'org.Hs.eg.db', # human gene annotation database
                 keyType = 'ENTREZID',
                 universe = names(y), # background gene set
                 ont = type,
                 readable = T) # return results with gene symbol annotation

  # Transform the results into a data frame
  return(go %>% 
           as.data.frame(.) %>% 
           dplyr::mutate(
             # Change gene ratio to numeric format
             GeneRatio = sapply(stringr::str_split(GeneRatio, "/"), 
                                function(y) as.numeric(y[1])/as.numeric(y[2])),
             # Change background ratio to numeric format
             BgRatio = sapply(stringr::str_split(BgRatio, "/"), 
                              function(y) as.numeric(y[1])/as.numeric(y[2])),
             Count = as.numeric(Count)))
}


# Function to create visualization for enrichment results
make_dotplot <- function(mydata, Count, type=c("GO","KEGG")){
  attach(mydata)
  # Get the size limits and breaks for the dot plot
  limits = c(min(Count),max(Count))
  modulus = mod(seq(limits[1], limits[2]), ceiling(length(seq(limits[1], limits[2])) / 10))
  
  breaks = if (max(Count) > 100) {
    unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 10, f= ceiling))
  } else if (max(Count) > 10 & max(Count) < 20){
    unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 2, f= ceiling))
  } else if (max(Count) > 20){
    unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 5, f= ceiling))
  } else {
    unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 1, f= ceiling))
  }
  detach(mydata)
  
  
  # Arrange the data frame based on the gene ratio
  if(type == "KEGG"){
    # get the top 10 KEGG pathways
    data <- mydata %>% 
      dplyr::arrange(desc(GeneRatio)) %>%
      dplyr::mutate(Description = fct_reorder(Description, GeneRatio)) %>%
      .[1:min(10, (dim(mydata)[1])),]
  } else {
    # get the top 10 GO terms for each category: BP, MF, CC
    data <- mydata %>%
      dplyr::arrange(desc(GeneRatio)) %>%
      dplyr::mutate(Description = fct_reorder(Description, GeneRatio)) %>%
      group_by(., ONTOLOGY) %>%
      dplyr::slice_head(n = 10)
  }
  
  # Create the dot plot
  plot <-         
    ggplot(data, aes(x = GeneRatio, y = Description, size = as.numeric(Count))) + 
    geom_point(aes(color = pvalue)) + 
    guides(color = guide_colorbar(title = "p-value", order = 1),
           size = guide_legend(title = "Gene count", order = 2)) +
    scale_size(range = c(2,7), limits = limits, breaks =  breaks) + 
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          strip.text.y = element_text(size = 14)) + 
    labs(size = "Count", y = "")
  
  # if it is a GO term, facet the plot based on the ontology
  if(type == "GO"){
    return(plot + 
             facet_grid(ONTOLOGY ~ ., scales = "free", 
                        labeller = as_labeller(
                          c("BP" = "Biological processes",
                            "MF" = "Molecular functions",
                            "CC" = "Cellular components"))) +
             theme(strip.text.y = element_text(size = 14))
    )
  } else {
    return(plot)
  }
}

# ---------------------- Semantic Similarity functions -------------------------
# Function to calculate semantic similarity between GO terms and return the score
# similarity matrix
make_simMatrix <- function(list){
  list(BP = calculateSimMatrix(list$BP$ID,
                               orgdb= org.Hs.eg.db,
                               ont="BP",
                               method="Wang"),
       CC = calculateSimMatrix(list$CC$ID,
                               orgdb= org.Hs.eg.db,
                               ont="CC",
                               method="Wang"),
       MF = calculateSimMatrix(list$MF$ID,
                               orgdb= org.Hs.eg.db,
                               ont="MF",
                               method="Wang"))
}

# Function to reduce GO terms with high collinearity based on their semantic 
# similarity scores
get_reducedTerms <- function(simm, scores, treshold){
  tmp <- reduceSimMatrix(simm,
                         scores,
                         threshold = treshold,
                         orgdb="org.Hs.eg.db")
  # Perform PCA on the reduced GO terms
  pca <- prcomp(simm)
  pca <- as.data.frame(pca$x[,1:2])
  pca <- pca %>% 
    setNames(c('x','y')) %>%
    rownames_to_column(., var = "go")
  
  tmp <- merge(tmp, pca, by = 'go')
  tmp$parentTerm <- as.factor(tmp$parentTerm)
  subset <- tmp [ tmp$termDispensability == 0, ]
  
  return(list(reduced = tmp, subset = subset))
}

# Function to visualise semantic similarity between GO terms on a PCA plot
make_GO_simplot <- function(reduced, subset){
  return(ggplot( data = reduced ) +
           geom_point( aes( x = x, y = y, colour = parentTerm, size = score), 
                       show.legend = F, alpha = I(0.6)) + 
           geom_point( aes( x = x, y = y, size = score), shape = 21, show.legend = F,
                       fill = "transparent", colour = I (alpha ("black", 0.6) )) +
           # Add confidence ellipses
           stat_ellipse(geom = "polygon", 
                        aes( x = x, y = y, colour = parentTerm, fill = parentTerm), 
                        type = "norm", level = 0.68, 
                        alpha = .3, linetype = 2, show.legend = F) + 
           scale_size( range=c(3, 10)) +
           theme_minimal() + 
           # Add labels for the representative GO terms
           geom_label_repel( data = subset, show.legend = F,
                             aes(x = x, y = y, label = term, colour = parentTerm), size = 3 ) +
           labs (y = "semantic space x", x = "semantic space y"))
}

make_GObase <- function (terms, genes) {
  # make sure that the gene symbols are comparable
  terms$geneID <- toupper(terms$geneID)
  genes$geneID <- toupper(genes$geneID)
  # create arrays of involved gene's names for each GO term
  tgenes <- strsplit(as.vector(terms$geneid), "/")
  # count the number of genes for each GO term
  count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
  # get the logFC values for each gene
  logFC <- sapply(unlist(tgenes), function(x) {
    genes$log2FoldChange[match(x, genes$geneID)]
  })
  # make sure that the logFC values are numeric
  if (class(logFC) == "factor") {
    logFC <- gsub(",", ".", gsub("\\.", "", logFC))
    logFC <- as.numeric(logFC)
  }
  s <- 1
  zsc <- c()
  for (c in 1:length(count)) {
    value <- 0
    e <- s + count[c] - 1
    value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 1, -1))
    # Z-score is the sum of the logFC values divided by the square root of the 
    # number of genes in the specific GO term
    zsc <- c(zsc, sum(value)/sqrt(count[c]))
    s <- e + 1
  }
  terms$zscore <- zsc
  return(terms)
}

# Create a complex plot of all GO categories and KEGG pathways and their activation
compound_GOplot <-function (data) {
  colnames(data) <- tolower(colnames(data))
  subset <- subset.data.frame(data, pvalue < 0.05)
  subset2 <- subset.data.frame(subset, slim == TRUE)
  # the dummy columns for background color
  dummy_col <- data.frame(category = factor(c("BP", "MF", "CC", "KEGG"),
                                            levels = c("BP", "MF", "CC", "KEGG")),
                          pvalue = data$pvalue[1:4], zscore = data$zscore[1:4],
                          size = 1:4,
                          generatio = 1:4)
  data <- data %>% 
    dplyr::mutate(category = factor(category, levels=c("BP", "MF", "CC", "KEGG")))
  
  (plot = ggplot(data) + 
      geom_point(aes(x = zscore, y = -log10(pvalue), size = generatio),
                 shape = 21, fill = "grey50", col = "black", alpha = .2) +
      geom_point(data = subset, 
                 aes(x = zscore, y = -log10(pvalue), fill = category, size = generatio),
                 shape = 21, col = "black", alpha = .5) +
      geom_rect(data = dummy_col, aes(fill = category),
                xmin = -Inf, xmax = Inf,
                ymin = -Inf, ymax = Inf,
                alpha = 0.1, show.legend = F) +
      geom_vline(xintercept = 0, lty = 2, col = "black") +
      facet_grid(cols = vars(category), #space = "free_x", scales = "free_x",
                 labeller = as_labeller(c("BP"="Biological processes",
                                          "MF"="Molecular functions",
                                          "CC"="Cell components",
                                          "KEGG"="KEGG pathways"))) +
      scale_size(range = c(5, 30), 
                 guide = guide_legend(
                   title = "Gene ratio",
                   position = "right",
                   override.aes = list(shape = 16, fill = "grey"),
                   order = 2)) + 
      scale_x_continuous(
        limits = c(-5, 5),
        expand = expansion(0.1), n.breaks = 10) + 
      labs(x = "z-score", 
           y = "-log (adj p-value)") + 
      scale_fill_manual(
        name = "Categories",
        labels = c("Biological processes","Molecular functions",
                   "Cell components","KEGG metabolic pathways"),
        values = c("BP" = "chartreuse4", 
                   "MF" = "brown2",
                   "CC" = "cornflowerblue",
                   "KEGG" = "orange"),
        aesthetics = c("fill","color"),
        guide = "none") + 
      geom_hline(yintercept = 1.3, col = "red") + 
      geom_label_repel(data = subset2, force = 0.5,
                       aes(x = zscore, y = -log10(pvalue), 
                           label = stringr::str_wrap(description, 35),
                           fill = category), alpha = .8,
                       fontface = 'bold',
                       box.padding = 5, max.overlaps = Inf,
                       arrow = arrow(length = unit(0.015, "npc")),
                       label.padding = .5, show.legend = F) +
      theme_bw(base_size = 14) + 
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 16),
            axis.line = element_line(colour = "grey80"),
            axis.ticks = element_line(colour = "grey80"),
            strip.text.x = element_text(size = 16),
            panel.border = element_rect(fill = "transparent", colour = "grey80"),
            plot.margin = margin(.5,1,.5,1, 'cm'),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = 14),
            legend.position = "bottom",
            panel.grid.major = element_line(colour = "grey80"),
            plot.background = element_rect(color = "white")))
}

# -------------------- Venn diagram & upset plots ------------------------------
# Create and save a simple Venn diagram
simple_venn <- function(set1, set2, set3, filename){
  # Extract the gene sets
  tmp <- list(
    set1[,c("geneID","log2FoldChange")],
    set2[,c("geneID","log2FoldChange")],
    set3[,c("geneID","log2FoldChange")]
  )
  names(tmp) <- names
  # Create a Venn diagram
  venn <- venn.diagram(x = list(tmp[[1]]$geneID, tmp[[2]]$geneID, tmp[[3]]$geneID),
                       category.names = names(tmp),
                       filename = file.path(date, plots_dir, paste0( filename, ".png")), 
                       output = T,
                       
                       # Output features
                       imagetype="png" ,
                       height = 960 , 
                       width = 960 , 
                       resolution = 300,
                       compression = "lzw",
                       # Circles
                       lwd = 2,
                       lty = 'blank',
                       fill = brewer.pal(3, "Pastel2"),
                       
                       # Numbers
                       cex = .6,
                       fontface = "bold",
                       fontfamily = "sans",
                       
                       # Set names
                       cat.cex = 0.6,
                       cat.fontface = "bold",
                       cat.default.pos = "outer",
                       cat.pos = c(-27, 27, 135),
                       cat.dist = c(0.055, 0.055, 0.085),
                       cat.fontfamily = "sans",
                       rotation = 1)
}

# Function to create a Venn diagram from the shared DEGs
make_vennbase <- function(set1, set2, set3, names){
  # Extract the gene sets
  tmp <- list(
    set1[,c("geneID","log2FoldChange")],
    set2[,c("geneID","log2FoldChange")],
    set3[,c("geneID","log2FoldChange")]
  )
  names(tmp) <- names
  # Extract the number of genes in each region of the venn diagram
  list <- GOVenn(tmp[[1]], tmp[[2]], tmp[[3]], plot = F)$table
  list <- lapply(list, rownames_to_column, var = "geneID")
  list[["A_only"]] <- list[["A_only"]] %>% dplyr::rename(logFC_A = logFC)
  list[["B_only"]] <- list[["B_only"]] %>% dplyr::rename(logFC_B = logFC)
  list[["C_only"]] <- list[["C_only"]] %>% dplyr::rename(logFC_C = logFC)
  
  # Create a data frame for the Venn diagram
  table <- do.call(rbind.fill,rev(list))
  table <- table %>% dplyr::distinct(geneID, .keep_all = T)
  
  return(list(venn = venn, table = table))
}

make_upsetbase <- function(table, vars){
  # Create an upset plot base: change log2FC values to logical (TRUE/FALSE)
  table <- table %>%
    dplyr::mutate(var1 = ifelse(is.na(logFC_A), FALSE, TRUE)) %>%
    dplyr::mutate(var2 = ifelse(is.na(logFC_B), FALSE, TRUE)) %>%
    dplyr::mutate(var3 = ifelse(is.na(logFC_C), FALSE, TRUE)) %>%
    dplyr::relocate(where(is.logical), .before = where(is.character))
  
  colnames(table)[1:3] <- vars
  
  return(table)
}

# ----------------------- Heatmap functions ------------------------------------
make_matrix <- function(mat, filter){
  # Calculate normalized expression values
  mat <- mat - rowMeans(mat)
  mat <- mat[which(row.names(mat) %in% filter),]
  
  return(mat)
  rm(mat)
}
