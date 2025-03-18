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
      dplyr::select(ONTOLOGY, ID, Description, Count, GeneRatio, pvalue, geneID) %>%
      dplyr::rename("category" = ONTOLOGY)
  })
  
  kegg <- kegg %>% 
    dplyr::select(ID, Description, Count, GeneRatio, pvalue, geneID) %>% 
    dplyr::mutate(category = "KEGG")
  df <- do.call(rbind, append(go, list("KEGG" = kegg)))
  df <- df %>% 
    dplyr::mutate(category = forcats::fct_relevel(category, "KEGG", after = Inf))
    
  return(df)
}

# A function to calculate TPM values based on a count matrix and feature lengths
get_TPM <- function(counts, effLen){
  n <- gsub("\\..*", "", names(effLen))
  g <- intersect(n, rownames(counts))
  
  effLen <- effLen[which(n %in% g)]
  x <- counts / effLen
  tpm.mat <- t( t(x) * 1e6 / colSums(x))
  
  return(tpm.mat)
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
  fdr <- fdrtool(tmp$stat, statistic = "normal", plot = F)
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
    dplyr::relocate(c("ensemblID","geneID","entrezID"), .after = everything()) #%>% 
    #dplyr::filter(complete.cases(.))
  
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
                     y = -log10(pvalue))) 
          + geom_point(mapping = aes(), color = "#c1c1c1", alpha = .3,
                       inherit.aes = T, size = 2.5, show.legend = F) 
          + geom_point(data = na.omit(sig), alpha = .7,
                       aes(color = significance), size = 2.5)
          + scale_color_manual(
            values = c(
              "Signif. up-regulated" = '#841f27',
              "Signif. down-regulated" = '#000f64'),
            guide = guide_legend(title = "Significance: ",
                                 override.aes = list(alpha = 1))
          )
          + labs(x = expression(paste(log[2], 'FoldChange')),
                 y = expression(paste(log[10], italic('P-value')))) 
          + theme(axis.title = element_text(size = 14), 
                  axis.text = element_text(size = 14), 
                  legend.position = 'none') 
          + scale_x_continuous(expand = expansion(0.2))
          # # Visualize log2FC threshold
          # + geom_vline(xintercept = c(-1.5, 1.5), 
          #              linetype = 'dotted', size = 1) 
          # # Visualize adjusted p-value threshold
          # + geom_hline(yintercept = -log10(0.05), 
          #              linetype = 'dotted', size = 1)
          # Add labels for significantly up- or down-regulated genes
          + geom_text(data = na.omit(sig), hjust = 0, vjust = 1.5, 
                      colour = 'black', position = 'identity', 
                      show.legend = F, check_overlap = T,
                      label = na.omit(sig)[,"geneID"]))
         + theme_minimal())
}

merge_vulcanplots <- function(p1, p2, p3, titles){
  legend <- get_legend(p1 +
                         guides(colour = guide_legend(title = "Significance: ")) + 
                         theme(legend.position="bottom", legend.box="horizontal",
                               legend.text = element_text(size = 12),
                               legend.title = element_text(size = 12)))
  
  return(plot <- cowplot::plot_grid(
    p1 + ggtitle(titles[1]) + theme(legend.position = "none"),
    p2 + ggtitle(titles[2]) + theme(legend.position = "none"),
    p3 + ggtitle(titles[3]) + theme(legend.position = "none"),
    ggplot() + theme_void(),
    legend,
    ggplot() + theme_void(),
    nrow = 2, ncol = 3, rel_widths = c(1,1,1), rel_heights = c(1,0.15)))
}

# ------------- Functions for overrepresentation analyses ----------------------
# Function to extract gene lists for KEGG pathway enrichment analysis
get_genelist <- function(.df, .filter){
  
  # Extract the background gene list of every expressed gene
  background <- .df %>%
    dplyr::mutate(effect_size=-log10(padj)*log2FoldChange) %>% 
    dplyr::arrange(desc(effect_size)) %>%
    dplyr::distinct(entrezID, .keep_all = T) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::pull("effect_size", name ="entrezID")
  # Extract the gene list of interest of DEGs
  interest <- .df %>%
    dplyr::mutate(effect_size=-log10(padj)*log2FoldChange) %>% 
    dplyr::filter(.filter) %>%
    dplyr::arrange(desc(effect_size)) %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::pull("effect_size", name ="entrezID")
  
  return(list(background = background, interest = interest))
}

get_zscore <- function(df, .expr, column = "geneID") {
  # make sure that the gene symbols are comparable
  genes <- toupper(df[[column]])
  .expr$geneID <- toupper(.expr$geneID)
  # create arrays of involved gene's names for each GO term
  tgenes <- strsplit(as.vector(genes), "/")
  # count the number of genes for each GO term
  count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
  # get the logFC values for each gene
  logFC <- sapply(unlist(tgenes), function(x) {
    if (!x %in% .expr$geneID) {
      return(0)
    } else {
      return(.expr$log2FoldChange[match(x, .expr$geneID)])
    }
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
    value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 1, ifelse( x < 0, -1, 0)))
    # Z-score is the sum of the logFC values divided by the square root of the 
    # number of genes in the specific GO term
    zsc <- c(zsc, sum(value)/sqrt(count[c]))
    s <- e + 1
  }
  df$zscore <- zsc
  return(df)
}


run_ora <- function(.interest, .background, .pathways){
  ora <- enricher(gene = names(.interest), # gene set of interest 
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  pAdjustMethod = "BH", 
                  universe = names(.background), # background gene set
                  TERM2GENE = dplyr::select(
                    .pathways,
                    gs_name,
                    human_entrez_gene))
  ora <- setReadable(ora, org.Hs.eg.db, keyType = "ENTREZID")
  return(list("ora" = ora))
}

extract_ora_results <- function(.ora, .db, .expr){
  .db <- .db %>% 
    dplyr::select(gs_name, gs_exact_source, gs_description) %>% 
    dplyr::distinct()
  # extract data frames
  df <- as.data.frame(.ora@result)
  # order on p-value
  df <- df[order(df$`p.adjust`, decreasing = F),]
  # change gene ratio to numeric values
  df <- df %>% 
    dplyr::select(!zScore) %>% 
    dplyr::mutate(GeneRatio = sapply(stringr::str_split(df$GeneRatio, "/"), 
                              function(x) 100*(as.numeric(x[1])/as.numeric(x[2]))),
           BgRatio = sapply(stringr::str_split(df$BgRatio, "/"), 
                            function(x) 100*(as.numeric(x[1])/as.numeric(x[2]))))
  # extract pathway IDs and description from database
  
  # add the pathway IDs and descriptions to the data frame
  ids <- df$ID
  database <- sapply(stringr::str_split(ids, "_"), 
                     function(x) return(x[1]))
  
  df$ID = .db[match(ids, .db$gs_name),][["gs_exact_source"]]
  df$Description = .db[match(ids, .db$gs_name),][["gs_description"]]
  df$Database = database
  
  # calculate activation z-score
  df = get_zscore(df, .expr)

  # extract significant results: adjusted p-value < 0.05
  sig_df <- df %>% 
    dplyr::filter(p.adjust < 0.05)
  
  #return data frames
  return(list("df" = df, "sig_df" = sig_df))
}

# III. GENE SET ENRICHMENT ANALYSIS
run_gsea <- function(.geneset, .terms){
  set.seed(42)
  ## run the GSEA analysis
  res <- GSEA(
    geneList = .geneset, # gene set of interest (ordered on effect size)
    minGSSize = 10, # minimum size of the gene set
    maxGSSize = 500, # maximum size of the gene set
    pvalueCutoff =1, # adjusted p-value cutoff
    eps =  0, # p-value cutoff (minimum)
    seed = TRUE, # seed for reproducibility
    pAdjustMethod = "BH", # p-value adjustment method
    TERM2GENE = dplyr::select(
      .terms,
      gs_name,
      human_entrez_gene
    ))
  res <- setReadable(res, org.Hs.eg.db, keyType = "ENTREZID")
  
  # extract data frame
  return(list("gsea" = res))
}

extract_gsea_results <- function(.gsea, .db, .expr){
  .db <- .db %>% 
    dplyr::select(gs_name, gs_exact_source, gs_description) %>% 
    dplyr::distinct()
  # extract data frames
  df <- as.data.frame(.gsea@result)
  # order on p-value
  df <- df[order(df$`p.adjust`, decreasing = F),]
  
  df <- df %>%
    dplyr::mutate(Direction = dplyr::case_when(NES > 0 ~ '+',
                                               NES < 0 ~ '-'))
  df <- df %>%
    dplyr::mutate(
      Database = sapply(stringr::str_split(ID, "_"), 
                        function(x) return(x[1])),
      Name = sapply(stringr::str_split(ID, "_"), 
                    function(x) return(paste(x[-1], collapse = "_"))))
  
  # add the pathway IDs and descriptions to the data frame
  ids <- df$ID
  df$ID = .db[match(ids, .db$gs_name),][["gs_exact_source"]]
  df$Description = .db[match(ids, .db$gs_name),][["gs_description"]]
  
  # calculate activation z-score
  df = get_zscore(df, .expr, column = "core_enrichment")
  
  # extract significant results: adjusted p-value < 0.05
  sig_df <- df %>% 
    dplyr::filter(p.adjust < 0.05)
  
  #return data frames
  return(list("df" = df, "sig_df" = sig_df))
}


cluster_enrichment <- function(.df, .db, column = "core_enrichment"){
  # 1.) Transform enriched terms to a long data table format
  table <- .df %>%
    dplyr::select(ID | tidyselect::starts_with(column)) %>%
    # Pivot the data frame to long format
    tidyr::pivot_longer(cols = tidyselect::starts_with(column),
                        values_to = "Genes", names_to = "Region") %>%
    tidyr::separate_rows(Genes, sep = "/") %>%
    dplyr::distinct(., .keep_all = T) %>% 
    # Group by geneID and concatenate the regions for genes in intersecting sets
    dplyr::mutate(Region = gsub(paste0(column, "."), "", Region)) %>%
    dplyr::group_by(ID, Genes) %>%
    dplyr::mutate(Region = paste(Region, collapse = "-")) %>% 
    # Filter out missing values
    dplyr::filter(!is.na(Genes)) %>% 
    dplyr::ungroup()
  
  # 2.) Calculate the similarity between enriched terms
  genes <- unique(table$Genes)
  total.pathway <- unique(table$ID)
  # Extract complete gene set of the background
  tmp <- list()
  for (p in total.pathway){
    tmp[[p]] <- .db %>%
      dplyr::filter(gs_exact_source == p) %>%
      dplyr::pull(gene_symbol)
  }
  # Create a kappa matrix 
  total.pathway <- tmp
  N <- length(total.pathway)
  kappa_mat <- matrix(0, nrow = N, ncol = N,
                      dimnames = list(names(total.pathway),
                                      names(total.pathway)))
  diag(kappa_mat) <- 1
  # Populate the kappa matrix 
  total <- length(genes)
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      genes_i <- total.pathway[[i]]
      genes_j <- total.pathway[[j]]
      both <- length(intersect(genes_i, genes_j))
      term_i <- length(base::setdiff(genes_i, genes_j))
      term_j <- length(base::setdiff(genes_j, genes_i))
      no_terms <- total - sum(both, term_i, term_j)
      observed <- (both + no_terms)/total
      chance <- (both + term_i) * (both + term_j)
      chance <- chance + (term_j + no_terms) * (term_i +
                                                  no_terms)
      chance <- chance/total^2
      kappa_mat[j, i] <- kappa_mat[i, j] <- (observed -
                                               chance)/(1 - chance)
    }
  }
  # Calculate similarity between pathways
  clu <- stats::hclust(stats::as.dist(1 - kappa_mat), method = "average")
  
  stats::heatmap(kappa_mat, distfun = function(x) stats::as.dist(1 - x), hclustfun = function(x) stats::hclust(x, method = "average"))
 
  kmax <- max(nrow(kappa_mat)%/%2, 2)
  if (kmax < 20) {
    kseq <- c(1:kmax)
  } else if (kmax < 50) {
    kseq <- c(2:19, seq(20, kmax %/% 10 * 10, 5))
  } else if (kmax < 200) {
    kseq <- c(2:19, seq(20, 49, 5), seq(50, kmax %/% 10 * 10, 10))
  } else {
    kseq <- c(2:19, seq(20, 99, 20), seq(100, kmax %/% 10 * 10, 50))
  }
  avg_sils <- c()
  avg_sils <- c()
  for (k in kseq) {
    avg_sils <- c(avg_sils, fpc::cluster.stats(stats::as.dist(1 - kappa_mat),
                                               stats::cutree(clu, k = k),
                                               silhouette = TRUE)$avg.silwidth)
  }
  k_opt <- kseq[which.max(avg_sils)]
  graphics::plot(clu)
  stats::rect.hclust(clu, k = k_opt)
  
  # Add clusters to the table
  clusters <- stats::cutree(clu, k = k_opt)
  clu_idx <- match(table[["ID"]], names(clusters))
  table$Cluster <- clusters[clu_idx]
  
  # 3.) Calculate the weight of each term
  linkage <- .df %>%
    #dplyr::inner_join(table[,c("ID","Cluster")], by = "ID") %>%
    dplyr::mutate(rank = 3 - rowSums(is.na(.)/3),
                  global.p = rowMeans(.df[,c(5,9,13)], na.rm = T)) %>% 
    dplyr::select(ID, rank, global.p) %>%
    dplyr::left_join(table, by = "ID") %>% 
    dplyr::group_by(Cluster, rank) %>% 
    dplyr::arrange(desc(rank), global.p) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(weight = rank*(-log10(global.p))*length(unlist(strsplit(Region, "-")))) 
  
  linkage <- linkage %>%
    dplyr::rowwise() %>%
    dplyr::mutate(weight = as.numeric(weight/max(linkage$weight))) %>%
    dplyr::ungroup() %>%
    dplyr::select(c("Genes", "ID", "weight", "Cluster")) %>%
    setNames(.,c("node1", "node2", "weight", "cluster")) %>%
    as.data.frame(.)
  
  
  net <- graph_from_data_frame(linkage)
  net <- simplify(net, remove.multiple = F, remove.loops = T)
  hs <- hub_score(net, scale = T, weights = linkage$weight)$vector
  
  
  linkage <- linkage %>%
    rowwise() %>% 
    dplyr::mutate(hub_score = hs[which(names(hs) == node1)])
  
  hub_scores <- linkage %>% 
    dplyr::group_by(node2, cluster) %>% 
    dplyr::summarise(hub_score = sum(hub_score, na.rm = T))
  
  .df = .df %>% 
    dplyr::right_join(., hub_scores, by = c("ID" = "node2"))

  representative.terms <- .df %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(hub_score)) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::pull(ID)

  .df <- .df %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(Type = ifelse(ID %in% representative.terms,
                                "Representative", "Member")) %>% 
    dplyr::ungroup(.)
  
  return(.df)
}

# Function to extract parent term of a GO term
get_parent <- function(x, y){
  parents <- GOMFPARENTS[[x]]
  if (length(parents) > 1){
    parents <- parents[which(parents %in% y)]
  }
  if (length(parents) == 0){
    parents <- x
  } # if the parent term is not in the list, return the term itself
  return(parents)
}

# Function to map GO terms to super families
map_to_super_family <- function(.df, .ref) {
  tmp <- .df
  i = 0
  c = 1
  while(c < 5 && i < 10){
    tmp <- tmp %>% 
      dplyr::group_by(GOID) %>%
      dplyr::mutate(count = n()) %>%
      dplyr::ungroup() %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(
        parent_GO = ifelse(count < 5,  get_parent(GOID, .ref), GOID)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(
        AnnotationDbi::select(GO.db, parent_GO, c("TERM"), "GOID")
      )
    i = i + 1
    print("Iteration: ", i)
    c <- min(tmp$count)
    print("Minimum term:", c)
  }
  table <- left_join(.df, tmp[,c("geneID","GOID","TERM")], by = "geneID",
                     suffix = c("",".parent")) %>% 
    dplyr::group_by(GOID.parent) %>%
    dplyr::mutate(count = n()) %>%
    dplyr::ungroup() %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(
      TERM.parent = ifelse(count < 5 | TERM.parent == "molecular_function",  "other", TERM.parent)) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(!count)
  
  table <- table %>% 
    mutate(TERM.parent = fct_relevel(TERM.parent, "other"))
  
  return(table)
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
  .ontology = names(list)
  .list = lapply(.ontology, function(x){
    calculateSimMatrix(list[[x]]$ID, orgdb =org.Hs.eg.db, ont=x, method="Wang")
  })
  .list = setNames(.list, names(list))
  return(.list)
}

# Function to reduce GO terms with high collinearity based on their semantic 
# similarity scores
get_reducedTerms <- function(simm, scores, limit){
  tmp <- reduceSimMatrix(simm,
                         scores,
                         threshold = 0.9,
                         orgdb="org.Hs.eg.db")
  # Perform PCA on the reduced GO terms
  pca <- prcomp(simm)
  pca <- as.data.frame(pca$x[,1:2])
  pca <- pca %>% 
    setNames(c('x','y')) %>%
    rownames_to_column(., var = "go")
  
  tmp <- merge(tmp, pca, by = 'go')
  tmp$parentTerm <- as.factor(tmp$parentTerm)

  subset <- tmp %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::arrange(desc(score)) %>%
    dplyr::slice_head(n = limit) %>% 
    dplyr::ungroup()
  
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

compare_GO <- function(.df1, .df2){
  data <- list()
  data <- lapply(list(.df1,.df2), function(x){
    return( x %>% 
              dplyr::select(!tidyselect::contains("enrichment")) %>%
              pivot_longer(cols = zscore.2D:Count.Tumor, 
                           names_pattern = "(.*)\\.(.*)$",names_to = c("measure", "condition"),
                           values_to = "value") %>%
              pivot_wider(names_from = measure, values_from = value) %>%
              dplyr::filter(complete.cases(.) & Type == "Representative") %>% 
              dplyr::mutate(condition = factor(condition, levels = c("2D", "3D","Tumor"),
                                               labels = c("Hypoxia (2D)", "Spheroid (3D)","Tumor"))) %>% 
              dplyr::group_by(condition) %>% 
              dplyr::arrange(desc(hub_score)) %>%
              dplyr::slice_head(., n = 7))
  })
  
  size_limits <- range(min(c(data[[1]]$Count, data[[2]]$Count)),
                       max(c(data[[1]]$Count, data[[2]]$Count)))
  fill_limits <- range(min(c(data[[1]]$zscore, data[[2]]$zscore)),
                       max(c(data[[1]]$zscore, data[[2]]$zscore)))
  
  p1 <- ggplot(data[[1]], 
               aes(x = condition, y = Description, fill = zscore, size = Count)) +
    geom_point(shape = 21, color = "black") + 
    scale_fill_gradient2(
        limits = fill_limits,
        low = "blue", mid = "white", high = "red", midpoint = 0,
        na.value = "grey50", guide = "colourbar", aesthetics = "fill"
      ) +
    scale_size(limits = size_limits,range = c(5,10)) + 
    guides(size = guide_legend(title = "Gene count", order = 1),
           fill = guide_colorbar(title = "Activation (z-score)", order = 2)) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) + 
    #facet_grid(Database~., drop = T, scales = "free_y") + 
    theme_bw() + 
    theme(legend.position = "right",
          legend.axis.line = element_line(colour = "#333333", size = 1.5),
          axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 16),
          axis.title = element_blank(),
          legend.title =  element_text(size = 16),
          legend.text = element_text(size = 16),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
  
  p2 <- ggplot(data[[2]], 
               aes(x = condition, y = Description, fill = zscore, size = Count)) +
    geom_point(shape = 21, color = "black") +
    scale_fill_gradient2(
        limits = fill_limits,
        low = "blue", mid = "white", high = "red", midpoint = 0,
        na.value = "grey50", guide = "colourbar", aesthetics = "fill"
      ) +
    scale_size(limits = size_limits,range = c(5,10)) +
    guides(size = guide_legend(title = "Gene count", order = 1),
           fill = guide_colorbar(title = "Activation (z-score)", order = 2)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 80)) +
    #facet_grid(Database~., drop = T, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "right",
          legend.axis.line = element_line(colour = "#333333", size = 1.5),
          axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 16),
          axis.title = element_blank(),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
  
  legend <- get_legend(p1)
  
  return(cowplot::plot_grid(
    p1 + theme(legend.position = "none"),
    p2 + theme(legend.position = "none"),
    legend,
    nrow = 1, ncol = 3, rel_widths = c(1,1,.3), rel_heights = c(1)))
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
simple_venn <- function(set1, set2, set3, names, filename){
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
  
  return(list(table = table))
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
# ----------------------- WGCNA analysis functions -----------------------------
# Function to perform WGCNA analysis
make_wgcna <- function(tpm, genes, power = 9, initial = TRUE){
  # Create a matrix for the WGCNA analysis
  matrix <- t(tpm[rownames(genes),])
  n_genes <- ceiling(ncol(matrix)/1000)*1000
  
  print(paste("Number of genes:", n_genes))
  
  if (initial){
    # Calculate the soft threshold
    powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
    
    soft = WGCNA::pickSoftThreshold(
      matrix,
      powerVector = powers,
      verbose = 5
    )
    
    df = data.frame(soft$fitIndices) %>%
      dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)
    
    plot = ggplot(df, aes(x = Power, y = model_fit, label = Power)) +
      # Plot the points
      geom_point() +
      # We'll put the Power labels slightly above the data points
      geom_text(nudge_y = 0.1) +
      # We will plot what WGCNA recommends as an R^2 cutoff
      geom_hline(yintercept = 0.80, col = "red") +
      # Just in case our values are low, we want to make sure we can still see the 0.80 level
      ylim(c(min(df$model_fit), 1.05)) +
      # We can add more sensible labels for our axis
      xlab("Soft Threshold (power)") +
      ylab("Scale Free Topology Model Fit, signed R^2") +
      ggtitle("Scale independence") +
      # This adds some nicer aesthetics to our plot
      theme_classic()
    
    return(list(vst = matrix, soft = df, plot = plot))
  } else {
    temp_cor <- cor       
    cor <- WGCNA::cor # Force it to use WGCNA cor function 
                      # (fix a namespace conflict issue)
    # Perform the WGCNA analysis
    wgcna <- WGCNA::blockwiseModules(matrix,
                                     # == Adjacency Function ==
                                     power = power,               
                                     networkType = "signed",
                                     # == Tree and Block Options ==
                                     deepSplit = 2,
                                     pamRespectsDendro = F,
                                     detectCutHeight = 0.75,
                                     minModuleSize = 30,
                                     maxBlockSize = n_genes,
                                     # == Module Adjustments ==
                                     reassignThreshold = 0,
                                     mergeCutHeight = 0.25,
                                     # == TOM == Archive the run results in TOM file (saves time)
                                     saveTOMs = T,
                                     saveTOMFileBase = "ER",
                                     # == Output Options
                                     numericLabels = T,
                                     verbose = 3)
    return(wgcna)
  }
  # Create a data frame for the WGCNA analysis
}

# Function to extract module genes and their expression
get_modules <- function(wgcna, plot_name){
  # Convert labels to colors for plotting
  mergedColors <- WGCNA::labels2colors(wgcna$colors)
  # Plot the dendrogram and the module colors underneath
  modules <- data.frame(
    gene_id = names(wgcna$colors),
    colors = WGCNA::labels2colors(wgcna$colors)
  )
  #
  svg(file.path(date, plots_dir, paste0(plot_name, ".svg")),
      width = 18, height = 10)
  WGCNA::plotDendroAndColors(
    wgcna$dendrograms[[1]],
    mergedColors[wgcna$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05)
  dev.off()
  
  return(modules)
}

module_expr <- function(expr, filter = c(''), colors){
  # Create a data frame for the module expression
  df <- expr %>% 
    dplyr::filter(! gene_id %in% filter)
  
  plot <- ggplot(df, aes(x = condition, y = averageExpr)) +
           geom_boxplot(aes(fill = module), position = "dodge") +
           geom_jitter(aes(group = module), color = "grey",
                       alpha = .5, position = position_jitter(width = 0.1)) +
           theme_bw() +
           scale_y_continuous(position = "right") +
           theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
           labs(title = "WGCNA Modules",
                x = "",
                y = "Average normalized expression")
  
  if (ncol(df) > 4){
    plot <- plot + 
      scale_fill_manual(
        name = "Cluster",
        values = colors) +
      facet_grid(rows = vars(df[,5]))
  } else {
    plot <- plot + 
      facet_wrap(~module)
  }
  return(plot)
}
# Function to visualize WGCNA module-trait correlations
module_corr <- function(module, colors, expression){
  # Get Module Eigengenes per cluster and reorder
  # them so similar modules are next to each other
  MEs0 <- WGCNA::orderMEs(module)
  
  
  plot <- draw(ComplexHeatmap::Heatmap(
    t(MEs0),
    name = "Module eigengenes",
    col = colorRamp2(c(min(MEs0), 0, max(MEs0)),
                     c("blue", "white", "red")),
    show_column_dend = T, show_column_names = T,
    show_row_names = T, show_row_dend = T,
    column_dend_reorder = F, row_split = 3,
    row_dend_reorder = F,
    clustering_distance_columns = "euclidean",
    column_title = "Module eigengenes",
    row_title = "Modules",
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.3f", t(dt)[i, j]), x, y, gp = gpar(fontsize = 5, angle = 90))},
    row_names_gp = gpar(fontsize = 8),
    column_title_gp = gpar(fontsize = 12),
    row_title_gp = gpar(fontsize = 12),
    bottom_annotation = heatmap.annot,
    heatmap_legend_param = list(title = "Module eigengenes"),
    width = unit(10, "cm"),
    height = unit(10, "cm")
  ))
  return(list(plot = plot))
}
# ----------------------- Heatmap functions ------------------------------------
make_matrix <- function(mat, filter){
  # Calculate normalized expression values
  mat <- mat - rowMeans(mat)
  mat <- mat[which(row.names(mat) %in% filter),]
  # mat <- mat[filter,]
  
  return(mat)
  rm(mat)
}

make_heatmap <- function(deseq, expr, dend, module, filter, coldata){
  # Create a expression matrix for the heatmap
  # matrix <- make_matrix(assay(deseq), 
  #                      module[filter, "gene_id"])
  # Create annotation for conditions and patients (columns)
  #column.annot <- HeatmapAnnotation(
  #  patient = as.factor(coldata[["patient"]]),
  #  condition = as.factor(coldata[["condition"]]),
  #  show_annotation_name = F,
  #  col =list(patient = c("593" = "navy", "673" = "khaki"),
  #            condition = c("normoxia.2D" = "steelblue",
  #                          "hypoxia.2D" = "salmon",
  #                          "physioxia.3D" = "red",
  #                          "physioxia.Tumour" = "darkred")))
  # Create annotation for WGCNA modules and module colors (rows)
  row.annot <- rowAnnotation(
    module = as.factor(module[filter, "colors"]),
    col = list(module = c(
      setNames(levels(as.factor(module[filter, "colors"])),
               levels(as.factor(module[filter, "colors"])))
    )),
    show_annotation_name = T, annotation_label = "Module colors", 
    width = unit(1, "in"))
  # Create the normalized expression heat map
  # hm1 <- Heatmap(matrix, name = "Normalized expression",
  #               cluster_rows = dend,
  #               show_column_dend = T, show_column_names = F,
  #               show_row_names = F, show_row_dend = F,
  #               cluster_columns = T, 
  #               clustering_distance_columns = "spearman",
  #               clustering_method_columns = "complete",
  #               column_dend_reorder = T,
  #               col = colorRamp2(c(min(matrix), 0, max(matrix)), 
  #                                c("green","black","red")),
  #               left_annotation = row.annot,
  #               bottom_annotation = column.annot)
  # Create the differential expression heat map
  plot <- Heatmap(as.matrix(expr), name = "log2FC",
                 na_col = "white", border = T,
                 cluster_rows = dend, column_dend_reorder = T,
                 show_column_dend = F, column_title = "Expression change",
                 show_column_names = T, show_row_names = F,
                 left_annotation = row.annot)
  # Merge the heat maps
  # plot <- draw(hm1 + hm2, merge_legend = T, padding = unit(c(15, 2, 2, 2), "mm"))
  
  return(list(plot = plot))
}



# ----------------------- SURFME functions -------------------------------------
# Function to perform SURFME analysis
get_CategoryExpressionPlot <- function(results, diff_genes, genes){
  subset <- subset.data.frame(results, subset = geneID %in% genes)
  colourPalette <- c('#c1c1c1', '#363636','#222222', '#000f64','#841f27')
  
  return(ggplot(data = na.omit(results)) 
         + geom_point(mapping = aes(x = log2FoldChange, y = -log10(padj), 
                                    colour = significance), size = 2.5) 
         + geom_point(data = subset, mapping = aes(x = log2FoldChange, y = -log10(padj)),
                      size = 2.5, fill = "transparent", colour = I (alpha ("yellow", 0.6) )) 
         + scale_color_manual(values = colourPalette) 
         + labs( x = expression(paste(log[2], 'FoldChange')),
                 y = expression(paste(log[10], italic('P')))) 
         + theme(axis.title = element_text(size = 14), 
                 axis.text = element_text(size = 14), 
                 legend.position = 'none') 
         + coord_cartesian(ylim = c(0, -log10(min(results$padj))))
         + geom_vline(xintercept = c(-(1.5), 1.5), linetype = 'dotted', size = 1) 
         + geom_hline(yintercept = -log10(0.05), linetype = 'dotted', size = 1) 
         + geom_label_repel(data = subset, aes(x = log2FoldChange, y = -log10(padj)),
                            colour = 'black', position = 'identity', max.overlaps = 20,
                            show.legend = F, label.padding = .5, direction = "both",
                            label = paste(subset$geneID)))
  
}