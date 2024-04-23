
# Recursively merges a list of data frames into a single data frame
merge.rec <- function(.list, ...){
  if(length(.list)==1) return(.list[[1]])
  Recall(c(list(merge(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
}


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
  df <- data.frame(tmp, geneID = rownames(assay(dds)), 
                   row.names = rownames(assay(dds))) %>% 
    dplyr::filter(complete.cases(.))
  
  # Add a column for significance using the specified thresholds
  significance <- rep('NS', nrow(df))
  significance[which(abs(df$log2FoldChange) > sig_log2FC & df$padj > sig_pval)] <- 'log2FoldChange'
  significance[which(abs(df$log2FoldChange) < sig_log2FC & df$padj < sig_pval)] <- '-Log10P'
  significance[which(df$log2FoldChange < -(sig_log2FC) & df$padj < sig_pval)] <- 'Signif. down-regulated'
  significance[which(df$log2FoldChange > sig_log2FC & df$padj < sig_pval)] <- 'Signif. up-regulated'
  df$significance <- as.factor(significance)
  # Create a data frame with only significant results
  sig_df <- df[which(df$significance == 'Signif. down-regulated' | 
                       df$significance == 'Signif. up-regulated'),]
  
  return(list(results = tmp, df = df, sig_df = sig_df))
  rm(list(tmp, fdr, df, significance, sig_df))
}


make_pca <- function(rld, group, labs, cols){
  tmp <- prcomp(t(assay(rld)), center =T, scale. = TRUE)
  plot <- ggbiplot(pcobj = tmp, choices = 1:2, scale = 1,
                   groups = rld@colData[[group]], var.axes = F, circle = T) +
    geom_point(size = 5, aes(color = rld@colData[[group]])) +
    geom_label_repel(aes(label = rld@colData@rownames, color = rld@colData[[group]]), 
                     size = 5, point.padding = 5, label.padding = .5, 
                     show.legend = F) + 
    scale_color_manual(name = "Conditions",
                       labels = labs,
                       values = cols) +
    theme(axis.text = element_text(size = 12, colour = "darkgrey"),
          axis.title = element_text(size = 12, face = "bold", colour = "black"),
          legend.text = element_text(size = 12, colour = "darkgrey"),
          legend.title = element_text(size = 12, face = "bold", colour = "black"))
  
  
  return(plot)
  #rm(list(tmp, plot))
}