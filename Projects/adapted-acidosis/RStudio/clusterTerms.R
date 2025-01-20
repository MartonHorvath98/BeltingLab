TOTAL.clusters <- do.call(rbind.fill, TOTAL.sets[c(5,6,7,9,13)])

enrichment.list <- list(U3017 = rbind(HGCC.GSEA$U3017$sig_df[,c("ID","core_enrichment")],
                                      HGCC.GO$U3017$sig_df[,c("ID","core_enrichment")]),
                        U3047 = rbind(HGCC.GSEA$U3047$sig_df[,c("ID","core_enrichment")],
                                      HGCC.GO$U3047$sig_df[,c("ID","core_enrichment")]),
                        U3054 = rbind(HGCC.GSEA$U3054$sig_df[,c("ID","core_enrichment")],
                                      HGCC.GO$U3054$sig_df[,c("ID","core_enrichment")]),
                        CCLD = rbind(CCLD.GSEA$sig_df[,c("ID","core_enrichment")],
                                     CCLD.GO$sig_df[,c("ID","core_enrichment")]))

TOTAL.clusters <- TOTAL.clusters %>% 
  dplyr::mutate(
    U3017.core_enrichment = left_join(., enrichment.list$U3017, by = "ID")$core_enrichment,
    U3047.core_enrichment = left_join(., enrichment.list$U3047, by = "ID")$core_enrichment,
    U3054.core_enrichment = left_join(., enrichment.list$U3054, by = "ID")$core_enrichment,
    CCLD.core_enrichment = left_join(., enrichment.list$CCLD, by = "ID")$core_enrichment) %>% 
  dplyr::select(c(1,2,3,4, # IDs
                  5,6,13, # U3017
                  7,8,14, # U3047
                  9,10,15, # U3054
                  11,12,16)) # CCLD

TOTAL.genes <- TOTAL.clusters %>%
  dplyr::select(!tidyselect::ends_with("FDR") & !tidyselect::ends_with("NES")) %>% 
  tidyr::pivot_longer(cols = tidyselect::ends_with("core_enrichment"),
                      values_to = "Genes", names_to = NULL) %>% 
  # Filter out missing values
  dplyr::filter(!is.na(Genes)) %>%
  # Pivot the data frame to long format, separate gene names
   tidyr::separate_rows(Genes, sep = "/") %>% 
  dplyr::distinct(., .keep_all = T)
  
  
total.genes <- unique(TOTAL.genes$Genes)
total.pathway <- unique(TOTAL.genes$ID)

tmp <- list()
for (pathway in total.pathway){
  tmp[[pathway]] <- rbind(terms,pathways) %>% 
    dplyr::filter(gs_exact_source == pathway) %>% 
    dplyr::pull(gene_symbol)
}

total.pathway <- lapply(tmp, unique)
N <- length(total.pathway)

kappa_mat <- matrix(0, nrow = N, ncol = N,
                    dimnames = list(names(total.pathway), 
                                    names(total.pathway)))
diag(kappa_mat) <- 1

total <- length(total.genes)
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

clu <- stats::hclust(stats::as.dist(1 - kappa_mat), method = "complete")

stats::heatmap(kappa_mat, distfun = function(x) stats::as.dist(1 - x), hclustfun = function(x) stats::hclust(x, method = "average")) 

#kmax <- max(nrow(kappa_mat)%/%2, 2)
kseq <- c(2:20, seq(50, 200, 50))
avg_sils <- c()

for (k in kseq) {
  avg_sils <- c(avg_sils, fpc::cluster.stats(stats::as.dist(1 - kappa_mat),
                                   stats::cutree(clu, k = k),
                                   silhouette = TRUE)$avg.silwidth)
}

plot(avg_sils)
k_opt <- kseq[which.max(avg_sils)]

graphics::plot(clu)
stats::rect.hclust(clu, k = 200)

clusters <- stats::cutree(clu, k = 200)
quantile(table(clusters))

clu_idx <- match(TOTAL.clusters[["ID"]], names(clusters))

TOTAL.clusters$Cluster <- clusters[clu_idx]

TOTAL.clusters <- TOTAL.clusters %>%
  dplyr::mutate(
    rank = 4 - rowSums(is.na(.)/3),
    
    global.p = rowMeans(TOTAL.clusters[,c(6,9,12,15)], na.rm = T)) %>% 
  dplyr::group_by(Cluster, rank) %>% 
  dplyr::arrange(desc(rank), global.p)


linkage <- dplyr::inner_join(TOTAL.clusters[,c("ID","Name","Cluster","rank","global.p")],
                             TOTAL.genes[,c("ID","Genes","Regions")], by = "ID") %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(weight = rank*(-log10(global.p))) 

linkage <- linkage %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(weight = as.numeric(weight/max(linkage$weight))) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(c("Genes", "Name", "weight", "Cluster")) %>% 
  setNames(.,c("node1", "node2", "weight", "cluster")) %>% 
  as.data.frame(.)

library(igraph)
library(viridis)
net <- graph_from_data_frame(linkage)
net <- simplify(net, remove.multiple = F, remove.loops = T) 

hs <- hub_score(net, scale = T, weights = linkage$weight)$vector

if(!"RCy3" %in% installed.packages()){
  BiocManager::install("RCy3", update = F)
}


TOTAL.genes <- TOTAL.genes %>% 
  dplyr::rowwise(.) %>% 
  dplyr::mutate(hub_score = hs[which(names(hs) == Genes)])

hub_scores <- TOTAL.genes %>% 
  dplyr::group_by(ID, Name) %>% 
  dplyr::summarise(core_enrichment = paste(Genes, collapse = "/"),
                   hub_score = sum(hub_score, na.rm = T))

TOTAL.clusters <- TOTAL.clusters %>% 
  dplyr::right_join(., hub_scores, by = c("ID","Name"))


representative.terms <- TOTAL.clusters %>%
  dplyr::group_by(Cluster) %>% 
  dplyr::arrange(desc(hub_score)) %>% 
  dplyr::slice_head(n = 1) %>% 
  dplyr::pull(Name)

TOTAL.clusters <- TOTAL.clusters %>% 
  dplyr::rowwise(.) %>% 
  dplyr::mutate(Type = ifelse(Name %in% representative.terms, 
                              "Representative", "Member"))

################################################################################
# 4. Vulcano visualization of the shared terms and pathways                    #
################################################################################
TOTAL.GSEA.vulcano <- list(
  U3017 = ungroup(TOTAL.clusters) %>% 
    dplyr::filter(grepl("U3017", Regions)) %>% 
    dplyr::select(ID, Name, U3017.NES, U3017.FDR, Cluster, Type) %>% 
    dplyr::rename(NES = U3017.NES, FDR = U3017.FDR),
  U3047 = ungroup(TOTAL.clusters) %>%
    dplyr::filter(grepl("U3047", Regions)) %>% 
    dplyr::select(ID, Name, U3047.NES, U3047.FDR, Cluster, Type)  %>% 
    dplyr::rename(NES = U3047.NES, FDR = U3047.FDR),
  U3054 = ungroup(TOTAL.clusters) %>%
    dplyr::filter(grepl("U3054", Regions)) %>% 
    dplyr::select(ID, Name, U3054.NES, U3054.FDR, Cluster, Type)  %>% 
    dplyr::rename(NES = U3054.NES, FDR = U3054.FDR),
  CCLD = ungroup(TOTAL.clusters) %>%
    dplyr::filter(grepl("CCLD", Regions)) %>% 
    dplyr::select(ID, Name, CCLD.NES, CCLD.FDR, Cluster, Type)  %>% 
    dplyr::rename(NES = CCLD.NES, FDR = CCLD.FDR))

interest_pathways <- read.csv("data/pathways-of-interest.txt", header = T,
                              sep = "\t", stringsAsFactors = T)

interest_clusters <- dplyr::left_join(interest_pathways, TOTAL.clusters, 
                                      by = c("ID", "Name")) %>% 
  dplyr::group_by(Cluster, Category) %>% 
  dplyr::select(Cluster, Category) %>%
  dplyr::distinct() 

cluster_palette = c("89"="#9ECAE1","63"="#4292C6","131"="#0D0887FF",
                    "146"="#4C02A1FF","186"="#7E03A8FF","2"="#A92395FF",
                    "90"="#FCBBA1","97"="#CC4678FF","167"="#E56B5DFF",
                    "159"="#EF3B2C","62"="#F89441FF","164"="#FDC328FF",
                    "147"="#F0F921FF")
                    
TOTAL.GSEA.vulcano.plots <- list()
TOTAL.GSEA.vulcano.plots <- lapply(TOTAL.GSEA.vulcano, function(x){
  plotClusters(.df = x, .pathways = interest_pathways, 
               .clusters = interest_clusters, .palette = cluster_palette)
})

TOTAL.GSEA.vulcano.plots <- sapply(names(TOTAL.GSEA.vulcano.plots), function(x){
  ggsave(file.path("Results",paste(x,"GSEA_Vulcano_plot",date,".png",sep = "_")),
         plot = TOTAL.GSEA.vulcano.plots[[x]], bg = "white",
         width = 20, height = 14, units = "in")
})

table <- TOTAL.GSEA.vulcano$U3017 %>% 
  dplyr::filter(Name %in% interest_pathways$Name) %>%
  tibble::column_to_rownames("Name") %>%
  dplyr::mutate(
    NES = round(NES, 4),
    FDR = round(FDR, 4)) %>%
  dplyr::mutate(FDR = case_when(
    FDR < 0.001 ~ paste("<0.001", "(***)"),
    FDR < 0.01 ~ paste(as.character(FDR), "(**)"),
    FDR < 0.05 ~ paste(as.character(FDR), "(*)"),
    TRUE ~ as.character(FDR),
  )) %>% 
  dplyr::select(!ID)

table_grob <- tableGrob(
  CCLD.enrichplot$table, 
  theme = ttheme_minimal(
    base_size = 14,
    core = list(
      fg_params = list(hjust = 0.5, x = 0.5, col = palette)),
    rowhead = list(
      fg_params = list(hjust = 0, x = 0,col = palette[c(5,1,2,3,4)]))
  )) 
