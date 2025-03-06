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
removeDuplicates <- function(.data, .column, .symbol){
  # Sort the data frame based on the log2FC values
  data <- .data[order(abs(.data[[.column]]), decreasing = TRUE),]
  # Remove duplicates
  data <- data[!duplicated(data[[.symbol]]),]
  return(data)
  
}

# III.) Differential expression analysis with Limma
limmaDEA <- function(.data, .design, .contrast){
  # Extract numeric columns
  mat <- .data[,sapply(.data, is.numeric)]
  colnames(mat) <- gsub("_\\(.*$", "\\1", colnames(mat))
  rownames(mat) <- .data$PROBEID
  # Create a design matrix
  t <- as.factor(.design)
  design <- model.matrix(~0 + t)
  colnames(design) <- levels(t)
  # Fit the linear model
  fit <- lmFit(mat, design)
  fit$genes$ID <- rownames(mat)
  fit$genes$Symbol <- .data$SYMBOL
  fit$genes$entrezID <- .data$ENTREZID
  
  # Fit the contrasts
  
  contrast.matrix <- makeContrasts(contrasts = .contrast, levels = t)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # Make deg table with getfitlimma
  res <- list()
  for (i in 1:ncol(contrast.matrix)){
    res[[i]] <- topTable(fit2, coef = i, number = Inf)
  }
  
  return(res)
}

# ---------------------------------------------------------------------- #
# Function - Illmina beadChip data preprocessing                         #
# ---------------------------------------------------------------------- #
# I.) Normalize transcript abundance using control probes
tryMapID <- function(ID, inp, outp){
  tryCatch(
    {
      suppressWarnings(
        mapIds(org.Hs.eg.db, keys = ID, column = outp, keytype = inp, multiVals = "first")
      )
    },
    error = function(e) {
      message(conditionMessage(e))
      return(NA)
    }
  )
}

switchAlias <- function(dbconn, ID){
  alias <- as.character(dbconn[which(dbconn$alias == ID),]$symbol)[1]
  
  if (is.na(alias)){
    return(ID)
  } else {
    return(alias)
  }  # returns TRUE
}


normalizeIllumina <- function(.df, .dbconn, .samples = colnames(.df$E)){
  x <- illuminaHumanv4ENTREZID
  
  #normalization and background correction
    df <- neqc(.df)
  # keep only selected samples
  df <- df[, .samples]
  #keep probes that are expressed in at least three arrays according to a detection p-values of 5%:
  expressed <- rowSums(df$other$Detection < 0.05) >= 3
  df <- df[expressed,]
  
  # get Illumina annotations
  annot <- df$genes
  annot <- cbind(PROBEID=rownames(annot), 
                 ENTREZID = sapply(as.list(x[rownames(annot)]), "[[", 1),
                 annot)
  annot <- annot %>%
    dplyr::select(!TargetID) %>%
    dplyr::filter(SYMBOL != "" & !is.na(SYMBOL)) %>% 
    dplyr::rowwise() %>% 
    # replace outdated gene symbols with the most recent ones
    dplyr::mutate(
      SYMBOL = ifelse(
        ENTREZID == "NA" | is.na(ENTREZID), no = SYMBOL, yes =  switchAlias(.dbconn, SYMBOL))
    )
  
  # map Illumina annotations against the org.Hs.eg.db database
  entrez <- tryMapID(annot$SYMBOL, "SYMBOL", "ENTREZID")
  # refine annotations missing from the Illumina database
  annot <- annot %>%
    rowwise() %>%
    dplyr::mutate(
      ENTREZID = ifelse(
        ENTREZID == "NA" | is.na(ENTREZID), no = ENTREZID, yes = entrez[[SYMBOL]]))
  
  # add annotation to the expression matrix
  df <- as.data.frame(df$E)
  df$PROBEID <- rownames(df)
  df <- merge(annot, df, by="PROBEID")
  # remove rows with missing ENTREZID
  idx <- which(is.na(df$ENTREZID))
  df <- df[-idx,]
  return(df)
}

# ---------------------------------------------------------------------- #
# Function - PCA                                                         #
# ---------------------------------------------------------------------- #
# I.) CREATE PCA PLOT
# The function to create a PCA plot with the ggbiplot package takes a prcomp
# matrix, groupings based on the metadata, the colors and names of the different
# treatments from the user input
plot_pca <- function(data, .groups, .labels = NULL, .values = NULL,
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
  return(plot)
}

# ---------------------------------------------------------------------- #
# Function - DIFF.EXP.                                                   #
# ---------------------------------------------------------------------- #

# I.) SET SIGNIFICANCE LEVELS
get_significance <- function(.df){
  return(.df %>% 
           # Rename data frame columns to make sense
           # setNames(., c("symbol", "log2FoldChange", "pvalue", "padj")) %>% 
           # Add a columns...
           dplyr::mutate(
             # ... for significance levels using the thresholds:
             #                            p-adj < 0.05, abs(log2FC) > 0.5
             significance = dplyr::case_when(
               abs(log2FoldChange) > 0.5 & padj > 0.05 ~ 'log2FoldChange',
               abs(log2FoldChange) < 0.5 & padj < 0.05 ~ 'Log10P',
               log2FoldChange < (-1)*0.5 & padj < 0.05 ~ 'Signif. down-regulated',
               log2FoldChange > 0.5 & padj < 0.05 ~ 'Signif. up-regulated',
               T ~ 'NS')) %>% 
           dplyr::filter(complete.cases(.)))
}

plot_vulcan <- function(.data, label = T){
  plot = ggplot(data = na.omit(.data), 
                aes(x = log2FoldChange, y = -log10(padj), colour = significance)) + 
    geom_point(mapping = aes(), inherit.aes = T, size = 2.5, alpha = 0.35) + 
    scale_color_manual(values = c("NS" = '#c1c1c1',
                                    "Log10P" = '#363636',
                                    "Log2FoldChange" = '#767676',
                                    "Signif. up-regulated" = '#841f27',
                                    "Signif. down-regulated" = '#000f64'),
                         name = "Significance") + 
    labs(x = expression(paste(log[2], 'FoldChange')),
           y = expression(paste(log[10], italic('FDR')))) +
    scale_x_continuous(expand = expansion(0.2)) + 
    # Visualize log2FC threshold
    geom_vline(xintercept = c(-0.5, 0.5), linetype = 'dotted', size = 1) +
    # Visualize adjusted p-value threshold
    geom_hline(yintercept = -log10(0.05), linetype = 'dotted', size = 1) +
    theme(axis.title = element_text(size = 14), 
            axis.text = element_text(size = 14), 
            legend.position = 'none')
  
  if (label == T){
    return(plot +
      # Add labels for significantly up- or down-regulated genes
      geom_text(data = subset(.data,
                                significance %in% c("Signif. up-regulated", 
                                                    "Signif. down-regulated")),
                  hjust = 0, vjust = 1.5, colour = 'black', position = 'identity', 
                  show.legend = F, check_overlap = T,
                  label = subset(.data,
                                 significance %in% c("Signif. up-regulated", 
                                                     "Signif. down-regulated"))[,"Symbol"]))
  } else {
    return(plot)
  }
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
      dplyr::select(c("Symbol","log2FoldChange")),
    .list[[2]] %>% 
      dplyr::filter(significance %in% c("Signif. up-regulated","Signif. down-regulated")) %>% 
      dplyr::select(c("Symbol","log2FoldChange")),
    .list[[3]] %>% 
      dplyr::filter(significance %in% c("Signif. up-regulated","Signif. down-regulated")) %>% 
      dplyr::select(c("Symbol","log2FoldChange"))
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
get_genelist <- function(.df, .filter, .value, .name){
  # Extract the background gene list of every expressed gene
  background <- .df %>%
    dplyr::distinct(entrezID, .keep_all = T) %>%
    dplyr::pull(.value, name = .name) %>% 
    sort(., decreasing = T)
  # Extract the gene list of interest of DEGs
  interest <- .df %>%
    dplyr::filter(.filter) %>%
    dplyr::distinct(entrezID, .keep_all = T) %>%
    dplyr::pull(.value, name = .name) %>% 
    sort(., decreasing = T)
  
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

# VI. PLOT GSEA RESULTS
# extract detailed enrichment statistics for a single term
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

# Plot running score for selected gene sets
plotRunningScore <- function(.df, .x, .y, .color, .palette){
  subset <- subset(.df, position == 1)
  
  return(ggplot(.df, aes(x = !!sym(.x)))
         + xlab(element_blank())
         + ylab("Enrichment score (ES)")
         + theme_bw(base_size = 14)
         + geom_hline(yintercept = 0, linetype = "dashed", color = "grey")
         + scale_x_continuous(expand = c(0, 0))
         + scale_y_continuous(expand = c(0.01, 0.01),
                              breaks = c(-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5))
         + geom_point(aes(x = !!sym(.x), y = !!sym(.y), color = !!sym(.color)),
                      show.legend = F, size = 1, data = subset)
         + scale_color_manual(values = .palette)
         + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.grid.major.y = element_blank(),
                 panel.grid.minor.y = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.y = element_text(colour = "black", size = 14))
  )
}

plotGeneRank <- function(.df, .x, .color, .facet, .palette){
  subset <- subset(.df, position == 1)
  facet <- as.formula(.facet)
  
  
  return(ggplot(.df, aes(x = !!sym(.x))) 
         + geom_linerange(aes(ymin = ymin, ymax = ymax, color = !!sym(.color)), 
                          show.legend = F, size = 1, data = subset) 
         + theme_bw(base_size = 14) 
         + xlab("Rank in ordered gene list") 
         + scale_x_continuous(expand = c(0, 0)) 
         + scale_y_continuous(expand = c(0, 0)) 
         + scale_color_manual(values = .palette) 
         + facet_grid(facet, scales = "free_y", switch = "y",
                      labeller = labeller(.default = Hmisc::capitalize)) 
         + theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(), 
                 panel.grid.major.y = element_blank(),
                 panel.grid.minor.y = element_blank(), 
                 axis.ticks.y = element_blank(), 
                 axis.text.y = element_blank(),
                 strip.placement = "outside",
                 strip.background = element_blank(),
                 strip.text.y.left = element_text(angle = 0, hjust = 0,
                                                  colour = "grey25", size = 14), 
                 axis.text.x = element_text(colour = "black", size = 14))
  )
}

getEnrichmentTable <- function(.df, .order, .name){
  return(.df %>%
           dplyr::arrange(.order) %>% 
           tibble::remove_rownames() %>% 
           tibble::column_to_rownames(.name) %>%
           dplyr::mutate(
             NES = round(NES, 4),
             FDR = round(p.adjust, 4)) %>%
           dplyr::mutate(FDR = case_when(
             FDR < 0.001 ~ paste("<0.001", "(***)"),
             FDR < 0.01 ~ paste(as.character(FDR), "(**)"),
             FDR < 0.05 ~ paste(as.character(FDR), "(*)"),
             TRUE ~ as.character(FDR),
           )) %>% 
           dplyr::select(c(NES, FDR)))
}



# ---------------------------------------------------------------------- #
# Function - Cluster analysis                                            #
# ---------------------------------------------------------------------- #
# a function calculating the Cohen's Kappa coefficient between a list of gene sets
cohen_kappa <- function(.list){
  N <- length(.list)
  kappa_mat <- matrix(0, nrow = N, ncol = N,
                      dimnames = list(names(.list), names(.list)))
  diag(kappa_mat) <- 1
  
  total <- length(unique(unlist(.list)))
  
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      genes_i <- .list[[i]]
      genes_j <- .list[[j]]
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
    }}
  return(kappa_mat)
}

# A function to create a network graph of gene sets based on the Cohen's Kappa coefficient
# then cluster the network based on interconnectivity using Louvain algorithm
get_cluster <- function(.df, .matrix, .column, .threshold = 0.25){
  # extract the list of enriched gene sets 
  genes <- split(strsplit(.df[["core_enrichment"]],"/"), .df[["ID"]])
  
  # subset the cohen's matrix for the gene sets of interest
  similarity.matrix <- .matrix[.df[['ID']], .df[['ID']]]
  
  # create an adjacency matrix based on the similarity matrix
  adjacency.matrix <- similarity.matrix > .threshold
  
  set.seed(42)
  # create a network graph from the adjacency matrix
  graph <- graph_from_adjacency_matrix(adjacency.matrix, mode = "undirected", diag = F)
  
  # annotate the nodes of the graph with the enrichment results
  V(graph)$Name <- .df$Name
  V(graph)$NES <- .df$NES
  V(graph)$geneRatio <- .df$geneRatio
  
  # Community detection (Louvain algorithm)
  clusters <- cluster_louvain(graph)
  # annotate nodes with cluster membership
  V(graph)$cluster <- membership(clusters)
  
  # extract cluster information
  cluster_summary <- data.frame(
    ID = names(V(graph)),
    cluster = as.factor(V(graph)$cluster))
  
  cluster_summary <- distinct(cluster_summary, .keep_all = T)
  
  df = .df %>% 
    dplyr::inner_join(., cluster_summary, by = "ID")
  
  return(list(graph = graph, df = df))
}

get_cluster_representative <- function(.cluster, .degs){
  # add gene expression values to the clusters
  linkage <- .cluster %>% 
    dplyr::select(ID, Name, core_enrichment, cluster) %>%
    tidyr::separate_rows(core_enrichment, sep = "/")
  
  linkage$logFC = .degs$log2FoldChange[match(linkage$core_enrichment, .degs$Symbol)]
  linkage$padj = .degs$padj[match(linkage$core_enrichment, .degs$Symbol)]
  
  # calculate normalized term weight
  linkage$weight = abs(linkage$logFC)*(-log10(linkage$padj))
  linkage$weight = linkage$weight/max(linkage$weight, na.rm = T)
  
  linkage = linkage %>% 
    dplyr::select(c("core_enrichment", "ID", "weight", "cluster")) %>%
    setNames(.,c("node1", "node2", "weight", "cluster")) %>% 
    as.data.frame(.)
  
  # create a network visualization of gene and GO-term relationships
  net <- graph_from_data_frame(linkage)
  net <- igraph::simplify(net, remove.multiple = F, remove.loops = T)
  
  # calculate hub score of each gene
  hs <-  igraph::hub_score(net, scale = T, weights = linkage$weight)$vector
  
  # summarize gene hub scores, to determine the final weight of each GO term
  linkage$geneRatio = .cluster$geneRatio[match(linkage$node2, .cluster$ID)]
  
  linkage <- linkage %>%
    dplyr::mutate(hub_score = geneRatio * hs[match(node1, names(hs))]) %>% 
    dplyr::rename(geneID = node1, ID = node2) %>% 
    dplyr::group_by(ID, cluster) %>%
    dplyr::summarise(core_enrichment = paste(geneID, collapse = "/"),
                     hub_score = sum(hub_score, na.rm = T))
  
  out_df <- inner_join(.cluster, linkage, by = c("ID","core_enrichment","cluster")) %>%
    dplyr::mutate(cluster = as.factor(cluster))
  
  # select cluster representatives according to calculated weight
  representative.terms <- out_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(hub_score)) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::pull(ID)
  
  out_df <- out_df %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(Representative = ifelse(ID %in% representative.terms, T, F))
  
  return(out_df)
}

filter_graph <- function(.graph, .threshold){
  cl <- table(V(.graph)$cluster)
  
  valid_cl <- which(cl >= .threshold)
  
  filtered_vertices <- which(V(.graph)$cluster %in% valid_cl)
  
  subgraph <- induced_subgraph(.graph, vids = filtered_vertices)
  
  return(subgraph)
}

plot_network <- function(.net, .layout, .df){
  edges = as.data.frame(as_edgelist(.net)) %>% 
    setNames(c("from", "to"))
  
  df = data.frame(x = as.data.frame(.layout)[,1],
                  y = as.data.frame(.layout)[,2],
                  ID = names(V(.net)),
                  Name = V(.net)$Description,
                  cluster = as.factor(V(.net)$cluster),
                  Representative = V(.net)$Representative,
                  weight = abs(.df[["hub_score"]][match(names(V(.net)), .df[["ID"]])]),
                  regulation = ifelse(.df[["NES"]][match(names(V(.net)), .df[["ID"]])] > 0, "UP", "DOWN")
  )
  
  lines = edges %>% 
    dplyr::mutate(
      from.x = df$x[match(edges$from, df$ID)],
      from.y = df$y[match(edges$from, df$ID)],
      to.x = df$x[match(edges$to, df$ID)],
      to.y = df$y[match(edges$to, df$ID)]
    )
  
  plot = ggplot() +
    stat_ellipse(data = df, geom = "polygon", 
                 aes(x = x, y = y, group = cluster), fill = "grey75", color = "grey15", #color = cluster, fill = cluster), 
                 type = "norm", level = 0.9,
                 alpha = .1, linetype = 2, show.legend = F) +
    geom_segment(data = lines, aes(x = from.x, y = from.y, xend = to.x, yend = to.y),
                 color = "grey25") +
    geom_point(data = df, aes(x = x, y = y, size = weight, fill = regulation), #fill = cluster, size = weight),
               shape = 21, colour = "black") +
    # geom_point(data = subset(df, Representative), aes(x = x, y = y), #fill = cluster, size = 5),
    #            shape = 21, colour = "orange", fill = "transparent", stroke = 2, size = 15) +
    geom_label_repel(data = df, 
                     aes(x = x, y = y,
                         label = ifelse(Representative, str_wrap(gsub("_"," ",Name), 30), "")),
                     max.overlaps = 100, min.segment.length = 0, size = 6,
                     color = "grey15", fill = "white", show.legend = F,
                     arrow = arrow(type = "closed", angle = 15, length = unit(0.1,"in")),
                     box.padding = unit(0.5,"in"), point.padding = unit(0.1,"in"),
                     force = 1, direction = "both") +
    scale_fill_manual(values = c("UP" = "darkred", "DOWN" = "blue"),
                      labels = c("UP" = "Activated", "DOWN" = "Inhibited"),
                      name = c("Regulation"),
                      guide = guide_legend(override.aes = c(size = 10))) +
    scale_size_continuous(guide = "none", range=c(5,15)) +
    theme_void() +
    theme(legend.title=element_text(size=28), 
          legend.text=element_text(size=28),
          legend.spacing=unit(1.5,"lines"),
          legend.position = "right")
  
  return(plot)
}

getCircplotData <- function(.cluster, .deg, .interest_cluster, .interest_cluster_genes, .palette){
  linkage <- .cluster %>% 
    dplyr::filter(ID %in% names(.interest_cluster)) %>% 
    dplyr::select(ID, Name, core_enrichment, cluster) %>%
    tidyr::separate_rows(core_enrichment, sep = "/") %>% 
    dplyr::left_join(.deg[,c("Symbol","log2FoldChange","padj")],
                     by = c("core_enrichment" = "Symbol")) %>% 
    dplyr::rowwise() %>%
    # ...as a function of the z-score and adjusted p-value
    dplyr::mutate(weight = abs(log2FoldChange)*(-log10(padj))) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(complete.cases(.)) %>% 
    dplyr::mutate(weight = as.numeric(weight/max(weight))) %>%
    dplyr::select(c("core_enrichment", "Name", "weight", "cluster")) %>%
    setNames(.,c("node1", "node2", "weight", "cluster")) %>%
    as.data.frame(.)
  
  # Transform input data in a adjacency matrix
  adjacencyData <- with(linkage, table(node1, node2))
  # circlize::chordDiagram(adjacencyData, transparency = 0.5)
  
  
  # create a network visualization of gene and GO-term relationships
  net <- graph_from_data_frame(linkage)
  net <- igraph::simplify(net, remove.multiple = F, remove.loops = T)
  # calculate hub score of each gene
  hs <-  igraph::hub_score(net, scale = T, weights = linkage$weight)$vector
  
  # summarize gene hub scores, to determine the final weight of each GO term
  linkage <- linkage %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(hub_score = hs[which(names(hs) == node1)])
  
  linkage_wide <- linkage %>% 
    dplyr::select(node1, node2, hub_score) %>%
    tidyr::pivot_wider(names_from = node2, values_from = hub_score) %>%
    dplyr::mutate(across(where(is.numeric), ~ifelse(is.na(.), 0, .))) %>% 
    tibble::column_to_rownames("node1") %>% 
    as.matrix()
  
  background = unique(linkage$node1[-c(which(linkage$node1 %in% .interest_cluster_genes))])
  focus = .interest_cluster_genes                   
  grid.col = c(.palette,
               setNames(rep("grey", length(background)), background),
               setNames(rep("red", length(focus)), focus))
  
  # Extract sector names
  from_sectors <- rownames(linkage_wide)
  to_sectors <- colnames(linkage_wide)
  
  border_mat = matrix(NA, nrow = nrow(linkage_wide), ncol = ncol(linkage_wide))
  border_mat[row.names(linkage_wide) %in% interest_cluster_genes, ] = "red"
  dimnames(border_mat) = dimnames(linkage_wide)
  
  return(list(data.mat = linkage_wide, border.mat = border_mat, grid.col = grid.col))
}

plotCircplot <- function(.path, .data, .color, .links, .labels){
  png(.path,
      width = 25, height = 15, res = 300, units = "in")
  circos.par(start.degree = 90)
  circlize::chordDiagram(.data, annotationTrack = "grid",
                         grid.col = .color, transparency = 0.5,
                         link.border = .links,
                         preAllocateTracks = list(track.height = 0.05))
  
  # Add labels only to the selected sectors
  circos.track(track.index = 1, panel.fun = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    
    # Label all 'to' sectors and only selected 'from' sectors
    if (sector_name %in% .labels) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    }
  }, bg.border = NA)
  dev.off()
  circos.clear()
}

# plot the enriched cluster on a vulcano plot of pathway's
plotClusters <- function(.df, .pathways){
  df <- .df %>% 
    dplyr::left_join(., .pathways, by = c("ID","Name")) %>%
    dplyr::mutate(Category = as.factor(Category))
  df <- df %>% 
    expand_grid(facet_var = unique(df$Category[!is.na(df$Category)]))
  
  return(
    ggplot() 
    + geom_point(data = df, aes(x = NES, y = -log10(p.adjust), fill = NES), 
                 size = 5, shape = 21, show.legend = F,
                 color = "black", alpha = 0.1) 
    # Visualize clusters of interest
    + geom_point(data = subset(df, !is.na(Category) & Category == facet_var),
                 aes(x = NES, y = -log10(p.adjust), colour = Category), 
                 show.legend = F, size = 5, alpha = 0.7)
    # Label the terms of interest
    + geom_label_repel(
      data = subset(df, ID %in% .pathways$ID & Category == facet_var),
      aes(x = NES, y = -log10(p.adjust), color = Category, 
          label = str_wrap(gsub("_"," ",Name), 25)),
      size = 5, max.overlaps = 100, min.segment.length = 0,
      arrow = arrow(type = "closed", angle = 15, length = unit(0.1,"in")),
      box.padding = unit(1.5,"in"), point.padding = unit(0.1,"in"),
      force = 1, direction = "both")
    + scale_color_manual(values = c("ECM"="#7E03A8FF",
                                    "HYPOXIA"="#F00000FF",
                                    "TGFb"="#09BFF9FF"))
    + scale_fill_gradient2(low = "blue", high = "red",
                           mid = "white", midpoint = 0)
    + facet_grid(~facet_var, scales = "free_x")
    + labs(x = expression('NES'),
           y = expression(paste(log[10], italic('FDR'))))
    + scale_x_continuous(expand = expansion(0.2))
    + theme_minimal() 
    + theme(axis.title = element_text(size = 14), 
            axis.text = element_text(size = 14, color = "#888888"),
            strip.text = element_blank(),
            # Controls spacing between facets
            panel.spacing = unit(5, "lines"),
            legend.position = 'none'))
}

score_plot <- function(.score, .formula, .ref_group, .x, .y, .title, number = F){
  # get formula
  formula <- as.formula(.formula)
  # statistical test
  stat.test = compare_means(formula, data = .score, 
                            method = "wilcox.test", p.adjust.method = "BH",
                            ref.group = .ref_group)
  
  range <- range(.score[[.y]])
  # Compute category counts
  counts <- .score %>%
    group_by(!!sym(.x)) %>%
    summarise(n = n()) %>%
    mutate(label = paste0(!!sym(.x), " (n=", n, ")"))  # Modify labels
  
  plot = ggplot(.score, aes(x=!!sym(.x), y=!!sym(.y))) + 
    #geom_jitter(aes(fill = Subtype), size = 5, alpha = 0.25) +
    geom_boxplot(color = "grey15", fill = "white",outlier.colour = "black",
                 outlier.shape = 21, outlier.stroke = 1.2, linewidth = .9,
                 alpha = 0.5) +
    scale_y_continuous(limits = range, expand = expansion(mult = c(0.1, 0.1)),
                       breaks = round(seq(floor(min(range)),
                                          ceiling(max(range)),
                                          2),1)) + 
    labs(y = str_wrap(.title, width = 25)) +
    theme_classic() + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "grey25"),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 14, color = "grey25")) + 
    stat_pvalue_manual(stat.test, hide.ns = TRUE, label = "p.signif",
                       size = 5,
                       remove.bracket = T, y.position = max(.score[[.y]]) - 0.1)
  
  if (number) {
    plot = plot +
      scale_x_discrete(labels = str_wrap(setNames(counts[["label"]], counts[[.x]]), 5))
  } else {
    plot = plot +
      scale_x_discrete(labels = str_wrap(counts[[.x]], 5))
  }
      
  return(list("plot" = plot, "stat.test" = stat.test))
}
