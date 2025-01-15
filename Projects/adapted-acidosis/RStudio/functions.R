# ---------------------------------------------------------------------- #
# Function - Microarray data processing                                  #
# ---------------------------------------------------------------------- #
# I.) Normailze transcript abundance
# The function normalizes the transcript abundance data using the RMA method
normalizeTranscript <- function(.data, .db){
  #rma does background normalization, log2 transformation and quantile normalization
  data <- rma(.data, target = "core")
  # annotating using info assembled by the Bioconductor Core Team
  data <- annotateEset(data, .db)
  
  # return expression matrix and add annotation
  df <- as.data.frame(exprs(data))
  df$PROBEID<- rownames(df)
  
  # Annotate gene names
  annot <- data@featureData@data
  df <- merge(df, annot, by="PROBEID")
  
  #Filter the gene expression to contain only values with genesymbol
  idx <- which(is.na(df$SYMBOL))
  df<- df[-idx,]
  
  return(df)
}

# II.) Remove duplicates
# Remove duplicate symbols based on lgFC
removeDuplicates <- function(.data, .column){
  # Sort the data frame based on the log2FC values
  data <- .data[order(abs(.data[[.column]]), decreasing = TRUE),]
  # Remove duplicates
  data <- data[!duplicated(data$SYMBOL),]
  return(data)

}

# III.) Differential expression analysis with Limma
limmaDEA <- function(.data, .design, .contrast){
  # Extract numeric columns
  mat <- .data[,sapply(.data, is.numeric)]
  colnames(mat) <- gsub("_\\(.*$", "\\1", colnames(mat))
  rownames(mat) <- .data$PROBEID
  print(head(mat))
  # Create a design matrix
  t <- as.factor(.design)
  design <- model.matrix(~0 + t)
  print(design)
  # Fit the linear model
  fit <- lmFit(mat, design)
  fit$genes$ID <- rownames(mat)
  fit$genes$Symbol <- .data$SYMBOL
  fit$genes$Entrez <- .data$ENTREZID
  
  # Fit the contrasts
  
  contrast <- makeContrasts(contrasts = .contrast, levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  # Make deg table with getfitlimma
  res <- decideTests(fit2, method="global", adjust.method="BH", p.value=0.05)
  degs <- topTable(fit2, number = "Inf")
  return(degs)
}

# IV.) Create DEG table
getFitLimma<-function (fit, results = NULL, digits = 3, adjust = "BH", method = "separate", 
                       F.adjust = "BH", scientificFormat=TRUE, sortByF=FALSE, filterByResults=FALSE,
                       ...) {
  if (!is(fit, "MArrayLM")) 
    stop("fit should be an MArrayLM object")
  if (!is.null(results) && !is(results, "TestResults")) 
    stop("results should be a TestResults object")
  if (is.null(fit$t) || is.null(fit$p.value)) 
    stop("results should include eBayes correction")
  method <- match.arg(method, c("separate", "global"))
  if(is.null(results))filterByResults=FALSE
  p.value <- as.matrix(fit$p.value)
  if (adjust == "none") {
    p.value.adj <- NULL
  } else {
    if (method == "separate"){
      p.value.adj <- p.value
      for (j in 1:ncol(p.value)) p.value.adj[, j] <- p.adjust(p.value[, j], method = adjust) 
    } else if (method == "global") {
      p.value.adj <- array(p.adjust(p.value, method = adjust),dim=dim(p.value))
      rownames(p.value.adj)<-rownames(p.value)
      colnames(p.value.adj)<-colnames(p.value)
    }
  }
  if (F.adjust == "none" || is.null(fit$F.p.value)){
    F.p.value.adj <- NULL 
  } else {
    F.p.value.adj <- p.adjust(fit$F.p.value, method = F.adjust)
  }
  rn <- function(x, digits = digits, fr=round){
    if (is.null(x)){
      NULL
    } else {
      if (is.matrix(x) && ncol(x) == 1) x <- x[, 1]
      fr(x, digits = digits)
    }
  }
  if(scientificFormat){
    fr <- function (x, digits = 3, ...) {
      x <- signif(x, digits)
      format(x, trim = TRUE, scientific = TRUE, ...)
    }
    dg1 <- digits
    dg2 <- digits
  } else {
    fr <- round
    dg1 <- digits + 2
    dg2 <- digits + 3
  }
  tab <- list()
  tab$Genes <- fit$genes
  tab$logFC <- rn(fit$coef, digits = digits)
  tab$Pv <- rn(p.value, digits = dg1, fr=fr)
  tab$PvAdj <- rn(p.value.adj, digits = dg2, fr=fr)
  # tab$F <- rn(fit$F, digits = digits)
  # tab$F.Pv <- rn(fit$F.p.value, digits = dg1, fr=fr)
  # tab$F.PvAdj <- tab$F.PvAdj <- rn(F.p.value.adj, digits=dg2, fr=fr)
  if(!is.null(results)){
    results<-unclass(results)
    junk<-sapply(colnames(results),function(nm){
      tab[[nm]]<<-results[,nm, drop=FALSE]
    })
  }
  tab <- data.frame(tab, check.names=FALSE, stringsAsFactors=FALSE)
  if(sortByF)tab<-tab[sort.list(fit$F.p.value),]
  colnames(tab)<-sub(" - ","-",colnames(tab))
  if(!is.null(colnames(fit$genes)))colnames(tab)[1:ncol(fit$genes)]<-colnames(fit$genes)
  if(filterByResults){
    signames<-rownames(results[rowSums(abs(results))>0, , drop=FALSE])
    tab<-tab[rownames(tab)%in%signames,]
  }
  return(tab)
}

# ---------------------------------------------------------------------- #
# Function - PCA                                                         #
# ---------------------------------------------------------------------- #
# I.) CREATE PCA PLOT
# The function to create a PCA plot with the ggbiplot package takes a prcomp
# matrix, groupings based on the metadata, the colors and names of the different
# treatments from the user input
plot_pca <- function(data, .groups, .labels, .values,
                     plot.ellipse = F, .ellipse){
  plot <- (
    # create a biplot object
    ggbiplot(pcobj = data,
             choices = 1:2, scale = 1, circle = T,
             groups = .groups, var.axes = F) +
      # add the dat apoints based on the first two PCAs
      geom_point(size = 5, aes(color = .groups, shape = .groups)) + 
      # define legend categories, names and colors
      scale_color_manual(name = "Groups",
                         labels = .labels,
                         values = .values) +
      # add a pretty theme
      theme(axis.text = element_text(size = 12, colour = "darkgrey"),
            axis.title = element_text(size = 12, face = "bold", colour = "black"),
            legend.text = element_text(size = 12, colour = "darkgrey"),
            legend.title = element_text(size = 12, face = "bold", colour = "black"),
            legend.position = "bottom", legend.direction = "horizontal"))
  
  if (plot.ellipse == T){
    plot <- plot + 
      # add confidence ellipses
      stat_ellipse(geom = "polygon", 
                     aes(color =  .ellipse, fill =  .ellipse), 
                     type = "norm", level = 0.68, alpha = .3, 
                     linetype = 2, show.legend = F)
  }
}

# ---------------------------------------------------------------------- #
# Function - DIFF.EXP.                                                   #
# ---------------------------------------------------------------------- #

# I.) SET SIGNIFICANCE LEVELS
get_significance <- function(.list){
  lapply(.list, function(x) {
    return(x %>% 
             # Rename data frame columns to make sense
             setNames(., c("symbol", "log2FoldChange", "pvalue", "padj")) %>% 
             # Add a columns...
             dplyr::mutate(
               # ... for ENTREZ gene identifiers, for downstream analyses
               entrezID = mapIds(org.Hs.eg.db, symbol, keytype = "SYMBOL", 
                                 column = "ENTREZID"),
               # ... for significance levels using the thresholds:
               #                            p-adj < 0.05, abs(log2FC) > 0.5
               significance = dplyr::case_when(
                 abs(log2FoldChange) > 0.5 & padj > 0.05 ~ 'log2FoldChange',
                 abs(log2FoldChange) < 0.5 & padj < 0.05 ~ 'Log10P',
                 log2FoldChange < (-1)*0.5 & padj < 0.05 ~ 'Signif. down-regulated',
                 log2FoldChange > 0.5 & padj < 0.05 ~ 'Signif. up-regulated',
                 T ~ 'NS')) %>% 
             dplyr::filter(complete.cases(.)))
  })
}

plot_vulcan <- function(.data){
  return(
    ggplot(data = na.omit(.data), 
         aes(x = log2FoldChange, y = -log10(padj), colour = significance))
    + geom_point(mapping = aes(), inherit.aes = T, size = 2.5, alpha = 0.35)
    + scale_color_manual(values = c("NS" = '#c1c1c1',
                                    "Log10P" = '#363636',
                                    "Log2FoldChange" = '#767676',
                                    "Signif. up-regulated" = '#841f27',
                                    "Signif. down-regulated" = '#000f64'),
                         name = "Significance")
    + labs(x = expression(paste(log[2], 'FoldChange')),
           y = expression(paste(log[10], italic('FDR'))))
    + scale_x_continuous(expand = expansion(0.2))
    # Visualize log2FC threshold
    + geom_vline(xintercept = c(-0.5, 0.5), linetype = 'dotted', size = 1)
    # Visualize adjusted p-value threshold
    + geom_hline(yintercept = -log10(0.05), linetype = 'dotted', size = 1)
    # Add labels for significantly up- or down-regulated genes
    + geom_text(data = subset(.data,
                              significance %in% c("Signif. up-regulated", 
                                                  "Signif. down-regulated")),
                hjust = 0, vjust = 1.5, colour = 'black', position = 'identity', 
                show.legend = F, check_overlap = T,
                label = subset(.data,
                               significance %in% c("Signif. up-regulated", 
                                                   "Signif. down-regulated"))[,"symbol"])
    + theme(axis.title = element_text(size = 14), 
            axis.text = element_text(size = 14), 
            legend.position = 'none'))
}

# ---------------------------------------------------------------------- #
# Function - SUBSETS                                                     #
# ---------------------------------------------------------------------- #
# I.) RECURSIVE MERGE FUNCTION
# Recursively merges a list of data frames into a single data frame
merge.rec <- function(.list, ...){
  if(length(.list)==1) return(.list[[1]])
  Recall(c(list(merge(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
}

# II.) EXTRACT GENE SETS
# The function collects genes that are representative to one or more conditions
# and returns in which region of the Venn diagram they will be found
get_regions <- function(.list, .names){
  # Extract the names and log2FC values of the significant genes
  tmp <- list(
    .list[[1]] %>% 
      dplyr::filter(significance %in% c("Signif. up-regulated","Signif. down-regulated")) %>% 
      dplyr::select(c("symbol","log2FoldChange")),
    .list[[2]] %>% 
      dplyr::filter(significance %in% c("Signif. up-regulated","Signif. down-regulated")) %>% 
      dplyr::select(c("symbol","log2FoldChange")),
    .list[[3]] %>% 
      dplyr::filter(significance %in% c("Signif. up-regulated","Signif. down-regulated")) %>% 
      dplyr::select(c("symbol","log2FoldChange"))
  )
  names(tmp) <- names(.list)
  
  # Extract the number of genes in each region of the Venn diagram
  sets <- GOVenn(tmp[[1]], tmp[[2]], tmp[[3]], plot = F)$table
  sets <- lapply(sets, rownames_to_column, var = "geneID")
  sets[["A_only"]] <- sets[["A_only"]] %>% dplyr::rename(logFC_A = logFC)
  sets[["B_only"]] <- sets[["B_only"]] %>% dplyr::rename(logFC_B = logFC)
  sets[["C_only"]] <- sets[["C_only"]] %>% dplyr::rename(logFC_C = logFC)
  
  # Create a data frame for the Venn diagram
  table <- do.call(rbind.fill,rev(sets))
  table <- table %>% dplyr::distinct(geneID, .keep_all = T)
  
  # Create an upset plot base: change log2FC values to logical (TRUE/FALSE)
  table <- table %>%
    dplyr::mutate("SET1" = ifelse(is.na(logFC_A), FALSE, TRUE)) %>%
    dplyr::mutate("SET2" = ifelse(is.na(logFC_B), FALSE, TRUE)) %>%
    dplyr::mutate("SET3" = ifelse(is.na(logFC_C), FALSE, TRUE)) %>%
    dplyr::relocate(where(is.logical), .before = where(is.character))
  colnames(table)[1:3] = .names
  
  # get the genes in each of the regions
  arranged <- arrange_venn(table, .names, extract_regions = T)
  arranged_df <- arrange_venn(table, .names, extract_regions = F)
  
  # return values: table and arranged sets
  return(list("table" = table, "sets" = arranged, "df" = arranged_df))
}

# III.) PLOT VENN DIAGRAM
plot_venn <- function(.data, .sets, .labels){
  intersection <- paste(.sets, collapse = "-")
  return((ggplot(.data)
          # --- BASE PLOT & LEGEND--- #
          # create the base plot with 1:1 aspect ratio and clear whit background
          + coord_fixed()
          + theme_void()
          # Add invisible gene sets to be able to add the legend to the plot
          + geom_point(aes(x = x, y = y, colour = Trend), 
               alpha = 0.01, size = 1)
          + scale_color_manual(
            labels = c("incoherent activation", "down-regulated", "up-regulated"),
            values = c("Change" = "orange", "DOWN" = "blue", "UP" = "darkred"),
            guide = guide_legend(override.aes = c(alpha = 1, size = 5)),
            aesthetics = "colour",
            ) 
          # --- SETS: BORDERS, FILL, LABELS --- #
          # Add the set borders of the Venn diagram
          + geom_venn_region(.data, sets =  .sets, 
                     alpha = 0.3, show.legend = F) 
          + geom_venn_circle(.data, sets =  .sets,
                     size = .5)
          # Fill every set white
          + scale_fill_venn_mix(.data, sets =  .sets,
                        guide='none', highlight = intersection,
                        inactive_color='NA', active_color = 'NA')
          # Add the names of the sets
          + geom_venn_label_set(.data, outwards_adjust = 2,
                        sets = .sets,
                        fill = alpha("black", .35), 
                        aes(label =  .labels))
          # --- DECORATE SETS; GENE EXPRESSION TRENDS --- #
          # add labels with numbers of incoherently activated genes
          + geom_venn_label_region(
            .data %>% dplyr::filter(Trend == "Change"),
            sets = .sets,
            aes(label=size ), fill = alpha("orange", .35),
            outwards_adjust=1.25, position=position_nudge(y=0.1))
          # add the numbers of up-regulated genes
          + geom_venn_label_region(
            .data %>% dplyr::filter(Trend == "UP"),
            sets = .sets,
            aes(label=size ), fill = alpha("darkred", .35),
            outwards_adjust=1.25, position=position_nudge(y=-0.1, x = 0.15))
          # add the numbers of down-regulated genes
          + geom_venn_label_region(
            .data %>% dplyr::filter(Trend == "DOWN"),
            sets = .sets,
            aes(label=size), fill = alpha("blue", .35),
            outwards_adjust=1.25, position=position_nudge(y=-0.1, x = -0.15)) 
          # --- SET THE THEME --- #
          + theme(
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 12))
  ))
  }

# ---------------------------------------------------------------------- #
# Function - PATHWAY ANALYSES                                            #
# ---------------------------------------------------------------------- #

# I. EXTRACT GENE LISTS
get_genelist <- function(.df, .filter){
  
  # Extract the background gene list of every expressed gene
  background <- .df %>%
    dplyr::mutate(effect_size=-log10(padj)*log2FoldChange) %>% 
    dplyr::arrange(desc(effect_size)) %>%
    dplyr::distinct(entrezID, .keep_all = T) %>%
    dplyr::pull("effect_size", name ="entrezID")
  # Extract the gene list of interest of DEGs
  interest <- .df %>%
    dplyr::mutate(effect_size=-log10(padj)*log2FoldChange) %>% 
    dplyr::filter(.filter) %>%
    dplyr::arrange(desc(effect_size)) %>%
    dplyr::distinct(entrezID, .keep_all = T) %>%
    dplyr::pull("effect_size", name ="entrezID")
  
  return(list(background = background, interest = interest))
}

# II. OVER-REPRESENTATION ANALYSIS
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

extract_ora_results <- function(.ora, .db){
  .db <- .db %>% 
    dplyr::select(gs_name, gs_exact_source, gs_description) %>% 
    dplyr::distinct()
  # extract data frames
  df <- as.data.frame(.ora@result)
  # order on p-value
  df <- df[order(df$`p.adjust`, decreasing = F),]
  # change gene ratio to numeric values
  df <- df %>% 
    mutate(GeneRatio = sapply(stringr::str_split(df$GeneRatio, "/"), 
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
    
  # extract significant results: adjusted p-value < 0.05
  sig_df <- df %>% 
    dplyr::filter(p.adjust < 0.1)
  
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
    pvalueCutoff = 1, # adjusted p-value cutoff
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




extract_gsea_results <- function(.gsea, .db){
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
  df$Database = .db[match(ids, .db$gs_name),][["gs_subcat"]]
  
  # extract significant results: adjusted p-value < 0.05
  sig_df <- df %>% 
    dplyr::filter(p.adjust < 0.1)
  
  #return data frames
  return(list("df" = df, "sig_df" = sig_df))
}

arrange_regions <- function(set1, set2, set3, set4, names){
  # GSEA on the GO terms from MSigDB
  tmp <- list(set1[,c("ID", "NES")],
              set2[,c("ID", "NES")],
              set3[,c("ID", "NES")],
              set4[,c("ID", "NES")])
  names(tmp) <- names
  # Merge the data frames
  table <- merge.rec(tmp, by = "ID", suffixes = c("",""), all = T) %>% 
    setNames(., c("ID", names))
  
  # Extract the allocation of each gene
  regions <- table %>%
    # Pivot the data frame to long format along the expression values
    tidyr::pivot_longer(cols = !ID, names_to = "Region", values_to = "Value") %>%
    # Filter out missing values
    dplyr::filter(!is.na(Value)) %>%
    # Group by geneID and concatenate the regions for genes in intersecting sets
    dplyr::group_by(ID) %>%
    dplyr::summarise(Regions = paste(Region, collapse = "")) %>%
    dplyr::ungroup() 
  
  # Merge the regions with the table
  table <- merge(table, regions, by = "geneID")
  
  # Calculate the number of genes in each subset
  regions <- regions %>% 
    dplyr::group_by(Regions) %>% 
    dplyr::summarise(Size = n()) %>% 
    dplyr::pull(Size, name = Regions)
  
  return(list(table = table, regions = regions)) 
}
