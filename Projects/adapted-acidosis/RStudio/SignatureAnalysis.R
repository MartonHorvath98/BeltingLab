################ IvyGap downloaded from GlioVis on 2024 05 24
# RNA-seq. The normalized count reads from the pre-processed data (sequence allignment and transcript abundance estimation) were log2 transformed after adding a 0.5 pseudocount (to avoid infinite value upon log transformation).
# Load Ivy Gap data
IvyGap <- read.table("./data/raw/IvyGap/2025-01-23_Ivy_GAP_expression.csv", header=TRUE,
                     sep = ",", dec = ".")
colnames(IvyGap) <- gsub(colnames(IvyGap), pattern = "X", replacement = "")
# Load Ivy Gap gene annotation
gene_IvyGap <- read.table("./data/raw/IvyGap/2025-01-23_Ivy_GAP_genes.csv", header=TRUE,
                          sep = ",", dec = ".")

IvyGap <- merge(gene_IvyGap[,c("gene_id", "gene_symbol")], IvyGap,
                by.x="gene_id", by.y="gene_id.rna_well_id")
rownames(IvyGap) <- IvyGap$gene_symbol
IvyGap <- IvyGap[,-c(1,2)]


# Load clinical data
clin1 <- read.table("./data/raw/IvyGap/2024-05-24_Ivy_GAP_pheno.txt", header=TRUE)
clin2 <- read.table("./data/raw/IvyGap/2025-01-23_Ivy_GAP_pheno.csv", header=TRUE,
                         sep = ",", dec = ".")
clin_IvyGap <- merge(clin1, clin2, by.x="rna_well_id", by.y="Sample")
rm(clin1,clin2)

# Get only primary tumors and samples that have location stated. Ivy gapd does not have IDH status information
# Out of 270 samples, 122 are primary GBM with Histology status 
Prim_IG <- clin_IvyGap[which(clin_IvyGap$Recurrence=="Primary" & !is.na(clin_IvyGap$Histology)),]
sel_IvyGap <- IvyGap[,which(colnames(IvyGap) %in% Prim_IG$rna_well_id)]

#Save RDatas
save(clin_IvyGap, IvyGap, Prim_IG, sel_IvyGap, file = "./RData/IvyGap_processedData.RData")

### Save data to excel
list_of_datasets <- list("PrimaryGBM_gexp" = sel_IvyGap, "PrimaryGBM_clinical" = Prim_IG,
                         "FullDataset_gexp" = IvyGap, "FullDataset_clinical" = clin_IvyGap)
write.xlsx(list_of_datasets, file = "./data/processed/IvyGap.xlsx", colNames=TRUE, rowNames=TRUE)

# ---------------------------------------------------------------------------- #
interest_genes <- read.csv("data/genes-of-interest.txt", header = T,
                          sep = "\t", stringsAsFactors = T)

Ivygap.GSVA <- list()
Ivygap.GSVA$gsva <- gsva(
  param = gsvaParam(as.matrix(sel_IvyGap),
                    list(fatlemon = as.character(interest_genes$Gene.name))),
  verbose = TRUE)
Ivygap.GSVA$zscores <- gsva(
  param = zscoreParam(as.matrix(sel_IvyGap),
                    list(fatlemon = as.character(interest_genes$Gene.name))),
  verbose = TRUE)
Ivygap.GSVA$df <- data.frame(
  sample = colnames(sel_IvyGap),
  gsva = as.numeric(Ivygap.GSVA$gsva),
  zscore = as.numeric(Ivygap.GSVA$zscores)
)

Ivygap.GSVA$df <- merge(Ivygap.GSVA$df, Prim_IG, by.x="sample", by.y="rna_well_id")
Ivygap.GSVA$df <- Ivygap.GSVA$df %>% 
  dplyr::rowwise(.) %>% 
  dplyr::mutate(
    across(Subtype.svm:Subtype_Verhaak_2010, ~gsub(x=.x,pattern="Proneural",replacement="Neural")),
    Subtype = paste(levels(factor(c_across(Subtype.svm:Subtype_Verhaak_2010))),
                    collapse = "/")
  ) %>% 
  dplyr::ungroup(.) %>% 
  dplyr::select(sample, gsva, zscore, Age, survival, status, Histology, 
                Subtype, MGMT_status, EGFR_amplification) %>% 
  dplyr::mutate(
    Histology = as.factor(Histology),
    Subtype = as.factor(Subtype),
    MGMT_status = as.factor(MGMT_status),
    EGFR_amplification = as.factor(EGFR_amplification)
  )


ggplot(Ivygap.GSVA$df) + 
  geom_boxplot(aes(x=Histology, y=gsva, fill=Histology)) +