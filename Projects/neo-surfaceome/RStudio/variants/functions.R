# ------------------- Data processing functions --------------------------------
# Recursively merges a list of data frames into a single data frame
merge.rec <- function(.list, ...){
  if(length(.list)==1) return(.list[[1]])
  Recall(c(list(merge(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
}
# A function to get the variants extracted and filtered through somatic variants
# and SURFME genes
merge_vcf <- function(.list, tool = "BCF"){
  data <- lapply(names(.list), function(x){
      return(.list[[x]][["fix"]] %>% 
               dplyr::select(c("CHROM","POS","REF","ALT","DP","QUAL")) %>% 
               dplyr::mutate(
                 # Create a unique identifier for each variant
                 variant_ID = stringr::str_glue("{CHROM}-{POS}-{REF}-{ALT}"),
                 # Determine the type of variant: SNP, DEL, INS, MNP
                 type = case_when(
                   nchar(REF) > nchar(ALT) ~ "DEL",
                   nchar(REF) < nchar(ALT) ~ "INS",
                   nchar(REF) == 1 & nchar(ALT) == 1 ~ "SNP",
                   T ~ "MNP"),
                 type = as.factor(type),
                 # Determine the class of SNPs: A>T, G>C, etc.
                 class = case_when(
                   type == "SNP" ~ stringr::str_glue("{REF}>{ALT}"),
                   T ~ NA_character_),
                 class = as.factor(class),
                 # Create an identifier for each sample
                 sample = as.factor(x)) %>% 
               dplyr::relocate(variant_ID, .before = "CHROM"))
  }) %>% 
      do.call(rbind, .)
  
  return(data)
}
  
summarise_vcf <- function(.df, conv, conv.class){
  
  data = .df %>% 
    dplyr::filter(!is.na(class))
  
  res = data %>%
    dplyr::mutate(class.con = conv[as.character(class)]) %>%
    dplyr::group_by(sample, method, class.con) %>% 
    dplyr::summarise(n = n(), .groups = "drop_last") %>%
    dplyr::mutate(fract = (n/sum(n)) * 100) %>%
    dplyr::mutate(TiTv = conv.class[class.con],
                  TiTv = factor(x = TiTv, 
                                levels = c("Tv", "Ti")))
  
  fract.classes = res %>% 
    reshape2::dcast(sample ~ class.con, mean, value.var = "fract")
  
  raw.classes = res %>% 
    reshape2::dcast(sample ~ class.con, mean, value.var = "n")
  
  titv = res %>%
    dplyr::group_by(sample, method, TiTv) %>%
    dplyr::summarise(n = sum(n), .groups = "drop") %>% 
    reshape2::dcast(sample ~ TiTv, mean, value.var = "n")
   
  return(list(res = res, fraction.contribution = fract.classes, raw.counts = raw.classes,
         TiTv.fractions = titv))
}
  

plot_summary <- function(df){
  # Plot 1 & Plot 3 uses the input data frame reordered by decreasing frequency
  df_1 = df %>% 
    dplyr::mutate(class.con = forcats::fct_reorder(class.con, desc(fract)))
  # Plot 2 summarizes the TiTv ratio for each sample
  df_2 = vcf.summary$res %>% 
    dplyr::group_by(sample, method, TiTv) %>% 
    dplyr::summarise(fract = sum(fract), n = sum(n)) %>% 
    dplyr::mutate(
      TiTv = forcats::fct_reorder(TiTv, desc(fract)))
  
  #return(list(df_1, df_2))
  # Pre-defined theme for pretty plots 
  
  my_theme <- theme(
      # plot background and grid
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "grey50", linetype = "longdash", size = 0.5),
      # axis and ticks
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.y = element_line(color = "black", size = 1),
      axis.ticks.y = element_line(color = "black", size = 1),
      axis.ticks.length.y = unit(0.3, "cm"), 
      # text elements
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      plot.margin = margin(.5, .5, .5, .5, "cm"))
  
  
  p1 <- ggplot(df_1, aes(x = class.con, y = fract)) +
    geom_point(aes(fill = class.con), color = "black", shape = 21, size = 3, alpha = .5, 
               position = position_jitter(width = 0.5),
               show.legend = F) +
    geom_boxplot(aes(fill = class.con), linewidth = 1, alpha = 0.5) + 
    labs(x = "", y = "% Mutations") + 
    scale_y_continuous(  # set y coordinate range and breaks
      limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
    scale_fill_manual(name = "Mutation type", values = mut_colors) +
    my_theme + 
    theme(legend.position = "none")
    
  
  p2 <- ggplot(df_2, aes(x = TiTv, y = fract)) +
    geom_point(color = "black", fill = "grey", shape = 21, size = 3, alpha = .5, position = position_jitter(width = 0.5)) +
    geom_boxplot(color = "black", fill = "grey25", linewidth = 1, alpha = 0.5) +
    labs(x = "", y = "") + 
    scale_y_continuous(  # set y coordinate range and breaks
      limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
    my_theme
  
  p3 <- ggplot(df_1, aes(x = sample, y = fract)) +
    #stacked bar of mutation classes
    geom_bar(stat = "identity", aes(fill = class.con), position = "stack") +
    facet_grid(.~method, scales = "free_y", labeller = labeller(
      method = c("bcftools" = "BCFtools", "gatk" = "GATK"))) +
    scale_fill_manual(name = "Mutation type",
                      values = mut_colors) +
    labs(x = "", y = "% Mutations") + 
    scale_y_continuous(  # set y coordinate range and breaks
      limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
    my_theme + 
    theme(
      panel.grid.major = element_blank(),
      # text elements
      axis.text.x = element_blank(),
      strip.text = element_text(size = 16),
      strip.background = element_rect(fill = "white")
    )
  
  legend <- get_legend(p3)
  
  # combine plots
  grid1 <- cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(2,1))#, align = "hv")
  grid2 <- cowplot::plot_grid(grid1, p3 + theme(legend.position = "none"),
                              nrow = 2, rel_heights = c(1.5, 1))#, align = "v")
  full.grid <- cowplot::plot_grid(grid2, legend, ncol = 2, rel_widths = c(1, 0.2))
  
  return(full.grid)
}

load_wgs <- function(path){
  tmp <-  read.csv(path, header = T, sep = "\t")
  
  split_info <- strsplit(tmp[["annotations"]], "\\|")
  # Convert the list to a matrix
  info_matrix <- do.call(rbind, lapply(split_info, "[", c(2:5, 18)))
  info_matrix <- data.frame(info_matrix)
  colnames(info_matrix) <- c("Consequence","Effect", "geneSymbol", "ensemblID", "rsID")
  
  df <- tmp %>% 
    tidyr::separate(., col = "variant_ID", into =c("CHROM", "POS", "REF", "ALT"),
                    sep = "-", remove = F) %>% 
    dplyr::select(c("CHROM", "POS", "REF", "ALT", "variant_ID"))
  
  df <- cbind(df, info_matrix)
  df <- df %>% 
    dplyr::mutate(across(where(is.character), as.factor))
  
  return(df)
}

filter_vcf <- function(df, somatic_table, surfme){
  som.df <- merge(df, somatic_table, by = c("CHROM", "POS", "REF","ALT", "variant_ID"))
  
  surf.df <- som.df[which(som.df$ensemblID %in% surfme),]
  
  return(list(df = df, somatic = som.df, surface = surf.df))
}

plot_top_genes <- function(.list, .ref, n, title){
  df <- do.call(rbind, lapply(.list, "[[", "somatic")) %>% 
    dplyr::group_by(geneSymbol) %>% 
    dplyr::summarise(variant_ID = paste(unique(variant_ID), collapse = "|")) %>% 
    dplyr::filter(geneSymbol != "") %>% 
    dplyr::mutate(n = stringr::str_count(variant_ID, "\\|") + 1) %>% 
    dplyr::arrange(desc(n)) %>% 
    dplyr::slice_head(n = n) %>% 
    tidyr::separate_rows(variant_ID, sep = "\\|") %>%
    dplyr::left_join(.ref[,c("variant_ID", "Consequence")],
                     by = c("variant_ID")) %>% 
    dplyr::mutate(geneSymbol = forcats::fct_reorder(geneSymbol, n)) %>% 
    dplyr::distinct(geneSymbol, variant_ID, Consequence, n)
  
  step <- if (max(df$n) >= 5) 5 else 2

  plot <- ggplot(data = df, aes(x = geneSymbol, fill = Consequence)) +
    geom_bar(stat = "count") + 
    scale_y_continuous(limits = c(0, max(df$n)), breaks = seq(0, max(df$n), max(df$n)%/%step), expand = c(0, 0)) +
    coord_flip() + 
    labs(title = title,
         x = "", y = "# Mutations") +
    scale_fill_manual(values = consequence_colors) +
    theme(
      # plot background and grid
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "grey50", linetype = "longdash", size = 0.5),
      # axis and ticks
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(color = "black", size = 1),
      axis.ticks.x = element_line(color = "black", size = 1),
      axis.ticks.length.x = unit(0.3, "cm"), 
      # text elements
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      plot.title = element_text(size = 20, hjust = 0.5, vjust = 2),
      plot.margin = margin(.5, .5, .5, .5, "cm")
    )
  return(list(df = df, plot = plot))
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