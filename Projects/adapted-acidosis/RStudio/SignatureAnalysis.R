lipid_droplet_genes <- read.csv("data/gene-signature-new.txt", header = T,
                           sep = "\t", stringsAsFactors = T)
subtype_hallmarks <- msigdbr_df %>% 
  dplyr::filter(gs_cat == "C2" & gs_subcat == "CGP" & 
                  grepl("VERHAAK_GLIOBLASTOMA_", gs_name)) %>%
  dplyr::pull(gs_name, name = gene_symbol) %>% 
  split(., sub("VERHAAK_GLIOBLASTOMA_", "", .))

subtype_hallmarks <- lapply(subtype_hallmarks, function(x) {
  return(names(x))
})
################ IvyGap downloaded from GlioVis on 2024 05 24
# RNA-seq. The normalized count reads from the pre-processed data (sequence allignment and transcript abundance estimation) were log2 transformed after adding a 0.5 pseudocount (to avoid infinite value upon log transformation).
# Load Ivy Gap data
IvyGap <- read.table("./data/raw/IvyGap/2025-01-23_Ivy_GAP_expression.csv", header=TRUE,
                     sep = ",", dec = ".")
IvyGap.genes <- read.table("./data/raw/IvyGap/2025-01-23_Ivy_GAP_genes.csv", header=TRUE,
                           sep = ",", dec = ".")
IvyGap <- merge(IvyGap.genes, IvyGap, by.x="gene_id", by.y="gene_id.rna_well_id")
row.names(IvyGap) <- IvyGap$gene_symbol
IvyGap <- IvyGap[,-c(1:5)]
colnames(IvyGap) <- gsub(x = colnames(IvyGap), pattern = "X", replacement = "")
# Load clinical data
clin1 <- read.table("./data/raw/IvyGap/2024-05-24_Ivy_GAP_pheno.txt", header=TRUE)
clin2 <- read.table("./data/raw/IvyGap/2025-01-23_Ivy_GAP_pheno.csv", header=TRUE,
                    sep = ",", dec = ".")
clin_IvyGap <- merge(clin2, clin1, by.x="rna_well_id", by.y="Sample")
rm(clin1,clin2)
# Get only primary tumors and samples that have location stated. Ivy gapd does not have IDH status information
# Out of 270 samples, 122 are primary GBM with Histology status
Prim_IG <- clin_IvyGap[which(clin_IvyGap$Recurrence=="Primary" & !is.na(clin_IvyGap$Histology)),]
sel_IvyGap <- IvyGap[,which(colnames(IvyGap) %in% Prim_IG$rna_well_id)]
Ivygap.score <- hack_sig(as.matrix(sel_IvyGap),
                         list(as.character(lipid_droplet_genes$SYMBOL)),
                         method = "zscore")
Ivygap.score <- setNames(as.data.frame(Ivygap.score), c("sample","zscore"))
Ivygap.score <- merge(Ivygap.score, Prim_IG, by.x="sample", by.y="rna_well_id")
Ivygap.score <- Ivygap.score %>%
  # dplyr::rowwise(.) %>%
  # dplyr::mutate(
  #   across(Subtype.svm:Subtype_Verhaak_2010, ~gsub(x=.x,pattern="Proneural",replacement="Neural"))
  # ) %>%
  # dplyr::ungroup(.) %>%
  dplyr::select(sample, zscore, Age, survival, status, Histology,
                Subtype.svm:Subtype_Verhaak_2010, MGMT_status, EGFR_amplification) %>%
  dplyr::mutate(
    Histology = factor(Histology,levels = c("Pseudopalisading cells", "Microvascular proliferation",
                                            "Leading Edge", "Infiltrating Tumor", "Cellular Tumor")),
    MGMT_status = as.factor(MGMT_status),
    EGFR_amplification = as.factor(EGFR_amplification)
  )
Ivygap.score$primary_subtype = apply(as.matrix(Ivygap.score[,7:10]), 1,function(x) unlist(Mode(as.factor(x)))[1])
Ivygap.type <- hack_sig(as.matrix(sel_IvyGap),
                        subtype_hallmarks,
                        method = "ssgsea", sample_norm = "separate", alpha = 0.5)
Ivygap.type <- Ivygap.type %>%
  tibble::column_to_rownames("sample_id") %>%
  dplyr::rename(Classical = CLASSICAL, Mesenchymal = MESENCHYMAL, Neural = NEURAL, Proneural = PRONEURAL)
Ivygap.type$Type = as.factor(colnames(Ivygap.type)[apply(Ivygap.type,1,which.max)])
Ivygap.score$gsea_subtype = Ivygap.type$Type
Ivygap.score <- Ivygap.score %>%
  dplyr::mutate(
    primary_subtype = factor(primary_subtype, levels = c("Mesenchymal","Classical","Proneural","Neural")),
    gsea_subtype = factor(gsea_subtype, levels = c("Mesenchymal","Classical","Proneural","Neural"))
  )
(Ivygap.plot$histology <- score_plot(.score = Ivygap.score, .formula = "zscore ~ Histology",
                                     .ref_group = "Pseudopalisading cells", .x = "Histology",
                                     .y = "zscore", .title = "Lipid droplet gene signature"))

# Save histology plot as png
ggsave(file.path("Results","signature", "IvyGap_signature_histology_log0_new_signature.png"), Ivygap.plot$histology$plot,
                 device = "png", width = 6, height = 4, units = "in")
# Save histology plot as svg
ggsave(file.path("Results", "signature","IvyGap_signature_histology_log0_new_signature.svg"), Ivygap.plot$histology$plot,
                 device = "svg", width = 6, height = 4, units = "in")

(Ivygap.plot$subtype <- score_plot(.score = Ivygap.score, .formula = "zscore ~ gsea_subtype",
                                  .ref_group = "Mesenchymal", .x = "gsea_subtype",
                                  .y = "zscore", .title = "Lipid droplet gene signature",
                                  number = T))

# Save subtype plot as png
ggsave(file.path("Results", "signature","IvyGap_signature_subtype_log0_new_signature.png"), Ivygap.plot$subtype$plot,
                 device = "png", width = 6, height = 4, units = "in")
# Save subtype plot as svg
ggsave(file.path("Results", "signature","IvyGap_signature_subtype_log0_new_signature.svg"), Ivygap.plot$subtype$plot,
                 device = "svg", width = 6, height = 4, units = "in")
write.xlsx(list("Histology" = Ivygap.plot$histology$stat.test,
                "Subtype" = CGGA.plot$subtype$stat.test),
           file.path("Results", "signature","IvyGap_histology+subtype_wilcox_test_new_signature.xlsx"),
           colNames=TRUE, rowNames=F)

rm(list = c("IvyGap", "clin_IvyGap", "Prim_IG", "sel_IvyGap", "Ivygap.type"))
# ---------------------------------------------------------------------------- #
################ IvyGap downloaded from GlioVis on 2025 02 20
# RNA-seq. The normalized count reads from the pre-processed data (sequence allignment and transcript abundance estimation) were log2 transformed after adding a 0.5 pseudocount (to avoid infinite value upon log transformation).
# Load Ivy Gap data
data_dir <- "./data/raw/TCGA_GBM/"
TCGA.expression <- read.table(file = file.path(data_dir, "TCGA_GBM_expression.txt"),
                              header = T, sep = "\t")
TCGA.clinical <- read.table(file = file.path(data_dir, "TCGA_GBM_pheno.txt"),
                            header = T, sep = "\t")

# exclude samples with missing gene expression data
row.names(TCGA.expression) <- TCGA.expression$Sample
TCGA.expression <- as.data.frame(t(TCGA.expression[,-1])) %>% 
  dplyr::mutate(across(everything(), ~ifelse(.x == -1, 0, .)))
TCGA.expression <- TCGA.expression[,apply(!is.na(TCGA.expression), 2, all)]

# exclude samples with missing gene expression data
TCGA.primary <- TCGA.clinical[which(TCGA.clinical$Sample %in% colnames(TCGA.expression)
                                    & TCGA.clinical$Recurrence=="Primary"  
                                    & TCGA.clinical$IDH1_status=="Wild-type"
                                    & !is.na(TCGA.clinical$Subtype)),]

sel_TCGA.expression <- TCGA.expression[,which(colnames(TCGA.expression) %in% TCGA.primary$Sample)]

#Save RDatas
save(TCGA.expression, TCGA.clinical, sel_TCGA.expression, TCGA.primary,
     file = "./RData/TCGA_processedData.RData")

### Save data to excel
list_of_datasets <- list("PrimaryGBM_gexp" = sel_TCGA.expression, "PrimaryGBM_clinical" = TCGA.primary,
                         "FullDataset_gexp" = TCGA.expression, "FullDataset_clinical" = TCGA.clinical)
write.xlsx(list_of_datasets, file = "./data/processed/TCGA.xlsx", colNames=TRUE, rowNames=TRUE)

TCGA.score <- hack_sig(as.matrix(sel_TCGA.expression),
                        list(as.character(lipid_droplet_genes$SYMBOL)), 
                        method = "zscore")

TCGA.score <- setNames(as.data.frame(TCGA.score), c("sample","zscore"))

TCGA.score <- merge(TCGA.score, TCGA.primary, by.x="sample", by.y="Sample")

TCGA.score <- TCGA.score %>% 
  dplyr::select(sample, zscore, Subtype_Verhaak_2010) %>% 
  dplyr::rename(primary_subtype = Subtype_Verhaak_2010)

TCGA.type <- hack_sig(as.matrix(sel_TCGA.expression),
                        subtype_hallmarks, 
                        method = "ssgsea", sample_norm = "separate", alpha = 0.5)

TCGA.type <- TCGA.type %>%
  tibble::column_to_rownames("sample_id") %>% 
  dplyr::rename(Classical = CLASSICAL, Mesenchymal = MESENCHYMAL, 
                Proneural = PRONEURAL, Neural = NEURAL)

TCGA.type$Type = as.factor(colnames(TCGA.type)[apply(TCGA.type,1,which.max)])
TCGA.score$gsea_subtype = TCGA.type$Type

TCGA.score <- TCGA.score %>% 
  dplyr::mutate(
    primary_subtype = factor(primary_subtype, levels = c("Mesenchymal","Classical","Proneural","Neural")),
    gsea_subtype = factor(gsea_subtype, levels = c("Mesenchymal","Classical","Proneural","Neural"))
  )

TCGA.plot <- list()
(TCGA.plot$subtype <- score_plot(.score = TCGA.score, .formula = "zscore ~ gsea_subtype",
                                  .ref_group = "Mesenchymal", .x = "gsea_subtype",
                                  .y = "zscore", .title = "Lipid droplet gene signature",
                                 number = T))

# Save subtype plot as png
ggsave(file.path("Results", "signature","TCGA_signature_subtype_log0_new_signature.png"), TCGA.plot$subtype$plot,
                 device = "png", width = 6, height = 4, units = "in")
# Save subtype plot as svg
ggsave(file.path("Results", "signature","TCGA_signature_subtype_log0_new_signature.svg"), TCGA.plot$subtype$plot,
                 device = "svg", width = 6, height = 4, units = "in")

write.xlsx(list("Subtype" = TCGA.plot$subtype$stat.test),
           file.path("Results", "signature","TCGA_subtype_wilcox_test_new_signature.xlsx"),
           colNames=TRUE, rowNames=F)

rm(list = c("TCGA.expression", "TCGA.clinical", "sel_TCGA.expression", 
            "TCGA.primary", "TCGA.type"))

# ---------------------------------------------------------------------------- #
################ CGGA downloaded from GlioVis on 2025 02 20
# RNA-seq. The normalized count reads from the pre-processed data (sequence allignment and transcript abundance estimation) were log2 transformed after adding a 0.5 pseudocount (to avoid infinite value upon log transformation).
# Load CGGA data
data_dir <- "./data/raw/CGGA/"
CGGA.expression <- read.table(file = file.path(data_dir, "CGGA_expression.txt"),
                              header = T, sep = "\t")
CGGA.clinical <- read.table(file = file.path(data_dir, "CGGA_pheno.txt"),
                            header = T, sep = "\t")

CGGA.primary <- CGGA.clinical[which(CGGA.clinical$Recurrence=="Primary"  
                                    & CGGA.clinical$IDH_status=="Wildtype"
                                    & CGGA.clinical$Histology=="GBM"),]

row.names(CGGA.expression) <- CGGA.expression$gene_name
CGGA.expression <- CGGA.expression[,-1]

CGGA.expression <- t(apply(CGGA.expression, 1, function(x) log2(x) + 1))
CGGA.expression <- as.data.frame(CGGA.expression) %>% 
  dplyr::select(which(colnames(CGGA.expression) %in% CGGA.primary$Sample)) %>% 
  dplyr::mutate(across(everything(), ~ifelse(is.infinite(.), 0, .)))

#Save RDatas
save(CGGA.expression, CGGA.clinical, CGGA.primary, file = "./RData/CGGA_processedData.RData")

### Save data to excel
list_of_datasets <- list("PrimaryGBM_gexp" = CGGA.expression, "PrimaryGBM_clinical" = CGGA.primary,
                         "FullDataset_gexp" = CGGA.expression, "FullDataset_clinical" = CGGA.clinical)

write.xlsx(list_of_datasets, file = "./data/processed/CGGA.xlsx", colNames=TRUE, rowNames=TRUE)

CGGA.score <- hack_sig(as.matrix(CGGA.expression),
                        list(as.character(lipid_droplet_genes$SYMBOL)), 
                        method = "zscore")

CGGA.score <- setNames(as.data.frame(CGGA.score), c("sample","zscore"))

CGGA.score <- merge(CGGA.score, CGGA.primary, by.x="sample", by.y="Sample")

CGGA.score <- CGGA.score %>% 
  dplyr::select(sample, zscore, Subtype) %>% 
  dplyr::rename(primary_subtype = Subtype)

CGGA.type <- hack_sig(as.matrix(CGGA.expression),
                        subtype_hallmarks, 
                        method = "ssgsea", sample_norm = "separate", alpha = 0.5)

CGGA.type <- CGGA.type %>%
  tibble::column_to_rownames("sample_id") %>% 
  dplyr::rename(Classical = CLASSICAL, Mesenchymal = MESENCHYMAL, Proneural = PRONEURAL, Neural = NEURAL)

CGGA.type$Type = as.factor(colnames(CGGA.type)[apply(CGGA.type,1,which.max)])
CGGA.score$gsea_subtype = CGGA.type$Type

CGGA.score <- CGGA.score %>% 
  dplyr::mutate(
    primary_subtype = factor(primary_subtype, levels = c("Mesenchymal","Classical","Proneural","Neural")),
    gsea_subtype = factor(gsea_subtype, levels = c("Mesenchymal","Classical","Proneural","Neural"))
  )

CGGA.plot <- list()
(CGGA.plot$subtype <- score_plot(.score = CGGA.score, .formula = "zscore ~ gsea_subtype",
                                  .ref_group = "Mesenchymal", .x = "gsea_subtype",
                                  .y = "zscore", .title = "Lipid droplet gene signature",
                                 number = T))

# Save subtype plot as png
ggsave(file.path("Results", "signature","CGGA_signature_subtype_log0_new_signature.png"), CGGA.plot$subtype$plot,
                 device = "png", width = 6, height = 4, units = "in")
# Save subtype plot as svg
ggsave(file.path("Results", "signature","CGGA_signature_subtype_log0_new_signature.svg"), CGGA.plot$subtype$plot,
                 device = "svg", width = 6, height = 4, units = "in")
write.xlsx(list("Subtype" = CGGA.plot$subtype$stat.test),
           file.path("Results", "signature","CGGA_subtype_wilcox_test_new_signature.xlsx"),
           colNames=TRUE, rowNames=F)
rm(list = c("CGGA.expression", "CGGA.clinical", "CGGA.primary", "CGGA.type"))


