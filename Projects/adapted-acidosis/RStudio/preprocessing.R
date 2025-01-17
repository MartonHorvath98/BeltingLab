# Created: 2024 05 21 ; Last Modified: 2025 01 16
# KGO & MH
################################################################################
#                              Data Processing                                 #
################################################################################
# Set working directory
wd <- getwd()
# Load packages
source(file = file.path(wd, "packages.R"))
# Source data processing functions
source(file = file.path(wd, "functions.R"))

# ---------------------------------------------------------------------------- #
# -   Svenja's CC +LD and CC no LD (Clariom D Human Pico) Affymetrix         - #
# ---------------------------------------------------------------------------- #
data_dir <- "./data/raw/ccLDvsccNoLD_Svenja_affy/"
CCLD.dat <- read.celfiles(list.celfiles(data_dir, full.names = TRUE))

# Transcript filtering
CCLD.expr <- normalizeTranscript(CCLD.dat, clariomdhumantranscriptcluster.db)
# log2Foldchange--> log2(exp)-log2(ctr)
CCLD.expr$log2FoldChange <- CCLD.expr$`CC_LD_(Clariom_D_Human).CEL` - CCLD.expr$`CC_noLD_(Clariom_D_Human).CEL`

####### Remove duplicate symbols based on lgFC
CCLD.expr <- removeDuplicates(.data = CCLD.expr,.column = "log2FoldChange",
                              .symbol = "SYMBOL")

####### Check how many up and down- regulated genes
length(which(CCLD.expr$log2FoldChange >=0.5)) # 6066
length(which(CCLD.expr$log2FoldChange <=(-0.5))) # 897

# save RData
dir.create("./RData", showWarnings = FALSE)
save(CCLD.expr, file = "./RData/CC_LD_noLD_processedData.RData")
# save file
dir.create("./data/processed", showWarnings = FALSE)
write.xlsx(CCLD.expr, file = "./data/processed/CC_LD_noLD_processedData.xlsx")

# ---------------------------------------------------------------------------- #
# -      Hugo's Primary cells 2D vs 3D (Clariom D Human Pico) Affymetrix     - #
# ---------------------------------------------------------------------------- #
data_dir <- "./data/raw/2Dvs3D_Hugo_primaryCells/"
HGCC.dat <- read.celfiles(list.celfiles(data_dir, full.names = TRUE))

# Transcript filtering
HGCC.expr <- normalizeTranscript(HGCC.dat, clariomdhumantranscriptcluster.db)

# Compare cell lines
HGCC.deg <- limmaDEA(.data = HGCC.expr,
                     .design = c("U3017_2D", "U3017_2D","U3017_2D",
                                 "U3017_3D", "U3017_3D", "U3017_3D",
                                 "U3047_2D", "U3047_2D","U3047_2D",
                                 "U3047_3D","U3047_3D","U3047_3D",
                                 "U3054_2D", "U3054_2D", "U3054_2D",
                                 "U3054_3D", "U3054_3D","U3054_3D"),
                     .contrast = c("tU3017_3D-tU3017_2D",
                                   "tU3047_3D-tU3047_2D",
                                   "tU3054_3D-tU3054_2D")) 
names(HGCC.deg) <- c("U3017", "U3047", "U3054")
# Compare global
HGCC.deg <- c(HGCC.deg,
              "global" = limmaDEA(.data = HGCC.expr,
                                  .design = c("2D","2D","2D","3D","3D","3D",
                                              "2D","2D","2D","3D","3D","3D",
                                              "2D","2D","2D","3D","3D","3D"),
                                  .contrast = c("t3D-t2D"))) 
# Remove duplicate IDs
HGCC.deg <- lapply(HGCC.deg, removeDuplicates, .column = "t", .symbol = "ID.Symbol")

HGCC.df <- merge.rec(HGCC.deg, by = c("ID.ID"), all = T, suffix = c("",""))
HGCC.df <- HGCC.df[,c(1:3, 5, # Identifier
                      4, 6:8, # U3017 2D vs 3D 
                      12, 14:16, # U3047 2D vs 3D 
                      20, 22:24, # U3054 2D vs 3D 
                      28, 30:32)] %>%  # Global 2D vs 3D
           setNames(., c("PROBEID", "Symbol", "EntrezID", "AveExpr",
             "log2FoldChange.U3017", "statT.U3017", "pvalue.U3017", "padj.U3017",
             "log2FoldChange.U3047", "statT.U3047", "pvalue.U3047", "padj.U3047",
             "log2FoldChange.U3054", "statT.U3054", "pvalue.U3054", "padj.U3054",
             "log2FoldChange.global", "statT.U3054", "pvalue.global", "padj.global"))
                  

####### Check how many up and down- regulated genes
length(which(HGCC.df$log2FoldChange.global >= 0.5)) # 2539
length(which(HGCC.df$log2FoldChange.global <= (-0.5))) # 1831

# Save RData
save(HGCC.deg, file = "./RData/HGCC_2D3D_processedData.RData")
# Save file
sapply(names(HGCC.deg), function(x){
  write.xlsx(HGCC.deg[[x]], file = paste0("./data/processed/HGCC_",x,"_2D3D_processedData.xlsx"))
})

# ---------------------------------------------------------------------------- #
# -     U87 Chronic Acidosis AA vs NA & HOX vs NOX (Illumina BeadChip)       - #
# ---------------------------------------------------------------------------- #
data_dir <- "./data/raw/U87_AAvsNA"
U87.dat <- read.ilmn(files=file.path(data_dir, "Sample_Probe_Summary.txt"),
                 ctrlfiles=file.path(data_dir, "Control_Probe_Summary.txt"))

samples <- c(
  # Chronic acidosis
  "200118400068_I", #13U87selctrlpH74 - Selection control pH 7.4
  "200118400035_I", #14U87selctrlpH74 - Selection control pH 7.4
  "200118400035_D", #15U87selctrlpH74 - Selection control pH 7.4
  "200118400035_K", #10U87selpH647 - Selection pH 6.4
  "200118400035_G", #11U87selpH647 - Selection pH 6.4
  "200118400035_F", #12U87selpH647 - Selection pH 6.4
  # Acute acidosis
  "200118400068_K", #7U87aactrl74 - Acute acidosis control pH 7.4
  "200118400068_D", #8U87aactrl74 - Acute acidosis control pH 7.4
  "200118400068_A", #9U87aactrl74 - Acute acidosis control pH 7.4
  "200118400068_C", #4U87aapH68 - Acute acidosis pH 6.8
  "200118400068_G", #5U87aapH68 - Acute acidosis pH 6.8
  "200118400033_E", #6U87aapH68 - Acute acidosis pH 6.8
  "200118400033_B", #1U87aapH64 - Acute acidosis pH 6.4
  "200118400068_B", #2U87aapH64 - Acute acidosis pH 6.4
  "200118400035_B", #3U87aapH64 - Acute acidosis pH 6.4
  # Normoxia
  "200118400033_G", #nox3 - atmospheric O2 (21%) control
  "200118400035_L", #nox2 - atmospheric O2 (21%) control
  "200118400033_H", #nox1 - atmospheric O2 (21%) control
  # Hypoxia
  "200118400035_E", #hox3 - 1% O2
  "200118400035_H", #hox2 - 1% O2
  "200118400033_I") #hox1 - 1% O2

#normalization and background correction
U87.expr <- normalizeIllumina(U87.dat, samples)

#Perform DEG analysis
U87.deg <- limmaDEA(.data = U87.expr,
                     .design = c("control_sel", "control_sel", "control_sel",
                                 "sel_pH647", "sel_pH647", "sel_pH647",
                                 "control_acu", "control_acu", "control_acu",
                                 "acu_pH68", "acu_pH68", "acu_pH68",
                                 "acu_pH64", "acu_pH64", "acu_pH64",
                                 "control_nox", "control_nox", "control_nox",
                                 "hypoxia", "hypoxia", "hypoxia"),
                     .contrast = c("tsel_pH647-tcontrol_sel",
                                   "tacu_pH68-tcontrol_acu",
                                   "tacu_pH64-tcontrol_acu",
                                   "thypoxia-tcontrol_nox"))

# Remove duplicate IDs
names(U87.deg) <- c("sel_pH647-control_sel", "acu_pH68-control_acu",
                    "acu_pH64-control_acu", "hypoxia-control_nox")
U87.deg <- lapply(U87.deg, removeDuplicates, .column = "t", .symbol = "ID.Symbol")

####### Check how many up and down- regulated genes
length(which(U87.deg$`sel_pH647-control_sel`$logFC >= 0.5)) # 1961
length(which(U87.deg$`sel_pH647-control_sel`$logFC <= (-0.5))) # 1945

# Save RData
save(U87.deg, file = "./RData/U87_AAvsNA_processedData.RData")
# Save file
sapply(names(U87.deg), function(x){
  write.xlsx(U87.deg[[x]], file = paste0("./data/processed/U87_",x,"_processedData.xlsx"))
})


#######    PANC1 AA vs NA (Clariom D Human Pico) Affymetrix
library(affycoretools)
library(oligo)
library(clariomdhumanprobeset.db)
library(clariomdhumantranscriptcluster.db)
library(pd.clariom.d.human)
library(openxlsx)
library(limma)

setwd("./RawData/PANC1_BEA21P006_affy/")
dat <- read.celfiles(list.celfiles())
setwd("E:/Lab/Collabs/Anna/Analysis2024")

######## Transcript filtering
#rma does background normalization, log2 transformation and quantile normalization
transcript.eset=rma(dat,target="core")

# annotating using info assembled by the Bioconductor Core Team
transcript.eset <- annotateEset(transcript.eset, clariomdhumantranscriptcluster.db)

#Get expression matrix and add annotation
PANC1_AAvsNA<- as.data.frame(exprs(transcript.eset))
PANC1_AAvsNA$PROBEID<- rownames(PANC1_AAvsNA)

annot_transc<- transcript.eset@featureData@data

PANC1_AAvsNA<- merge(PANC1_AAvsNA, annot_transc, by="PROBEID")

#Filter the gene expression to contain only values with genesymbol
idx<- which(is.na(PANC1_AAvsNA$SYMBOL))
PANC1_AAvsNA<- PANC1_AAvsNA[-idx,]


###### DEA
# Do differential expression to all samples
gexp<- PANC1_AAvsNA[,1:7]
rownames(gexp)<- PANC1_AAvsNA$PROBEID
gexp<-gexp[,-1]

t <- factor(c("PANC1_NA", "PANC1_NA","PANC1_NA",
              "PANC1_AA", "PANC1_AA", "PANC1_AA"))
design <- model.matrix(~0+t)

library(limma)
fit <- lmFit(gexp,design)
fit$genes$ID <- rownames(gexp)
fit$genes$Symbol <- PANC1_AAvsNA$SYMBOL
fit$genes$Entrez <- PANC1_AAvsNA$ENTREZID

#--- Contrasts of control groups against tested groups
contrasts <- makeContrasts(tPANC1_AA - tPANC1_NA,
                           levels=design) 

ct.fit <- eBayes(contrasts.fit(fit, contrasts))
res.fit<-decideTests(ct.fit,method="global", adjust.method="BH", p.value=0.05)

#--- Make deg table with getfitlimma
source("./Scripts/getFitLimma.R")
deg.limma<-getFitLimma(ct.fit,res.fit)

PANC1_AAvsNA_deg<-merge(PANC1_AAvsNA[,c(1,8,10,9,2:7)], deg.limma[,c(1,4:6,10)], by.x="PROBEID", by.y="Genes.ID")

PANC1_AAvsNA_deg<-PANC1_AAvsNA_deg[,c(1:4,11:14,5:10)]

rm(list = setdiff(ls(), c("PANC1_AAvsNA_deg")))

#remove duplicates based on the probe with highest abs FC
symbols<- as.character(unique(PANC1_AAvsNA_deg[which(duplicated(PANC1_AAvsNA_deg$SYMBOL)),"SYMBOL"]))

result<- NULL
for (i in 1:length(symbols)) {
  gene<- symbols[i]
  mat<- PANC1_AAvsNA_deg[which(PANC1_AAvsNA_deg$SYMBOL==gene),]
  mat<-mat[order(abs(mat$logFC), decreasing = TRUE),]
  result<- rbind(result, mat[1,])
}

PANC1_AAvsNA_deg<- PANC1_AAvsNA_deg[-which(PANC1_AAvsNA_deg$SYMBOL%in%symbols),]
PANC1_AAvsNA_deg<- rbind(PANC1_AAvsNA_deg, result)

####### Check how many up and down- regulated genes
length(which(PANC1_AAvsNA_deg$logFC>=0.5)) # 1040
length(which(PANC1_AAvsNA_deg$logFC<=(-0.5))) # 661

#Save RDatas
save(PANC1_AAvsNA_deg, file = "./RData/PANC1_AAvsNA_processedData.RData")
#write.xlsx(PANC1_AAvsNA_deg, file = "./ProcessedData/PANC1_AAvsNA_processedData.xlsx")





################ IvyGap downloaded from GlioVis on 2024 05 24
# RNA-seq. The normalized count reads from the pre-processed data (sequence allignment and transcript abundance estimation) were log2 transformed after adding a 0.5 pseudocount (to avoid infinite value upon log transformation).
# Load Ivy Gap data
IvyGap<- read.table("./RawData/IvyGap/2024-05-24_Ivy_GAP_expression.txt", header=TRUE)
# Load clinical data
clin_IvyGap<-read.table("./RawData/IvyGap/2024-05-24_Ivy_GAP_pheno.txt", header=TRUE)

# Get only primary tumors and samples that have location stated. Ivy gapd does not have IDH status information
# Out of 270 samples, 122 are primary GBM with Histology status 
Prim_IG<- clin_IvyGap[which(clin_IvyGap$Recurrence=="Primary" & !is.na(clin_IvyGap$Histology)),]
sel_IvyGap<-merge(Prim_IG, IvyGap, by="Sample")

#Save RDatas
save(clin_IvyGap, IvyGap, Prim_IG, sel_IvyGap, file = "./ProcessedData/IvyGap_processedData.RData")


# Save as xlsx
Samples<-IvyGap$Sample

IvyGap<-IvyGap[,-1]
IvyGap<-t(IvyGap)
IvyGap<-as.data.frame(IvyGap)
colnames(IvyGap)<-Samples


Samples<-sel_IvyGap$Sample

sel_IvyGap<-sel_IvyGap[,c(18:ncol(sel_IvyGap))]
sel_IvyGap<-t(sel_IvyGap)
sel_IvyGap<-as.data.frame(sel_IvyGap)
colnames(sel_IvyGap)<-Samples


list_of_datasets <- list("PrimaryGBM_gexp" = sel_IvyGap, "PrimaryGBM_clinical" = Prim_IG,
                         "FullDataset_gexp" = IvyGap, "FullDataset_clinical" = clin_IvyGap)
write.xlsx(list_of_datasets, file = "./ProcessedData/IvyGap.xlsx", colNames=TRUE, rowNames=TRUE)


