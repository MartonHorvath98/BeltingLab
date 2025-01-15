# Created: 2024 05 21 ; Last Modified: 2025 01 15
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
CCLD.expr <- removeDuplicates(CCLD.expr, "log2FoldChange")

medians <- rowMedians(as.matrix(CCLD.expr[,c(2,3)]))
hist(medians, 100, col = "cornsilk1", freq = FALSE, 
     main = "Histogram of the median intensities", 
     border = "antiquewhite4",
     xlab = "Median intensities")
####### Check how many up and down- regulated genes
length(which(CCLD.expr$log2FoldChange >=0.5)) # 6082
length(which(CCLD.expr$log2FoldChange <=(-0.5))) # 867

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
HGCC.res <- limmaDEA(.data = HGCC.expr,
                     .design = c("U3017_2D", "U3017_2D","U3017_2D",
                                 "U3017_3D", "U3017_3D", "U3017_3D",
                                 "U3047_2D", "U3047_2D","U3047_2D",
                                 "U3047_3D","U3047_3D","U3047_3D",
                                 "U3054_2D", "U3054_2D", "U3054_2D",
                                 "U3054_3D", "U3054_3D","U3054_3D"),
                     .contrast = c("tU3017_3D-tU3017_2D",
                                   "tU3047_3D-tU3047_2D",
                                   "tU3054_3D-tU3054_2D")) 
# Compare global
HGCC.global <- limmaDEA(.data = HGCC.expr,
                        .design = c("2D", "2D","2D",
                                    "3D", "3D", "3D",
                                    "2D", "2D","2D",
                                    "3D","3D","3D",
                                    "2D", "2D", "2D",
                                    "3D", "3D","3D"),
                        .contrast = c("t3D-t2D")) 

HGCC.deg <- dplyr::inner_join(HGCC.expr, HGCC.res,
                             by = c("PROBEID" = "Genes.ID",
                                    "ENTREZID" = "Genes.Entrez",
                                    "SYMBOL" = "Genes.Symbol")) %>% 
  dplyr::inner_join(., HGCC.global,
                    by = c("PROBEID" = "Genes.ID",
                           "ENTREZID" = "Genes.Entrez",
                           "SYMBOL" = "Genes.Symbol")) %>% 
  dplyr::mutate(across(starts_with("Pv"), as.numeric)) %>%
  dplyr::select(c(1, 20:22, # Identifier
                  23, 26, 29, # U3017 results
                  24, 27, 30, # U3047 results
                  25, 28, 31, # U3054 results
                  35:37, # Global results
                  2:19)) # Expression values

colnames(HGCC.deg)[5:16] <- paste(c("log2FoldChange","p.value","padj"), 
                                 rep(c(".U3017", ".U3047", ".U3054", ""), each = 3),
                                 sep = "")
HGCC.deg <- removeDuplicates(HGCC.deg, "log2FoldChange")

rm(list = c("HGCC.res", "HGCC.global"))

####### Check how many up and down- regulated genes
length(which(HGCC.deg$log2FoldChange >= 0.5)) # 2543
length(which(HGCC.deg$log2FoldChange <= (-0.5))) # 1836

# Save RData
save(HGCC.deg, file = "./RData/HGCC_2D3D_processedData.RData")
# Save file
write.xlsx(HGCC.deg, file = "./data/processed/HGCC_2D3D_processedData.xlsx")

###############################################################################
##      U87 AA vs NA (Illumina BEadChip)
####### Chronic Acidosis   
library(illuminaHumanv4.db)
library(limma)
library(ggplot2)
library(reshape2)
library(viridis)

setwd("./RawData/U87_AAvsNA/")
x <- read.ilmn(files="./Sample_Probe_Summary.txt",ctrlfiles="./Control_Probe_Summary.txt")
setwd("E:/Lab/Collabs/Anna/Analysis2024")

#Wanted barcodes - Chronic Acidosis
#13U87selctrlpH74 - 200118400068_I - Selection control pH 7.4
#14U87selctrlpH74 - 200118400035_I - Selection control pH 7.4
#15U87selctrlpH74 - 200118400035_D - Selection control pH 7.4
#10U87selpH647 - 200118400035_K - Selection pH 6.4
#11U87selpH647 - 200118400035_G - Selection pH 6.4
#12U87selpH647 - 200118400035_F - Selection pH 6.4

##Wanted barcodes - Acute Acidosis
#7U87aactrl74 - 200118400068_K - Acute acidosis control pH 7.4
#8U87aactrl74 - 200118400068_D - Acute acidosis control pH 7.4
#9U87aactrl74 - 200118400068_A - Acute acidosis control pH 7.4

#4U87aapH68 - 200118400068_C - Acute acidosis pH 6.8
#5U87aapH68 - 200118400068_G - Acute acidosis pH 6.8
#6U87aapH68 - 200118400033_E - Acute acidosis pH 6.8

#1U87aapH64 - 200118400033_B - Acute acidosis pH 6.4
#2U87aapH64 - 200118400068_B - Acute acidosis pH 6.4
#3U87aapH64 - 200118400035_B - Acute acidosis pH 6.4

samples<- c("200118400068_I", "200118400035_I", "200118400035_D",
            "200118400035_K", "200118400035_G", "200118400035_F",
            "200118400068_K", "200118400068_D", "200118400068_A",
            "200118400068_C", "200118400068_G", "200118400033_E",
            "200118400033_B", "200118400068_B", "200118400035_B")


#normalization and background correction
y <- neqc(x)

#keep probes that are expressed in at least three arrays according to a detection p-values of 5%:
expressed <- rowSums(y$other$Detection < 0.05) >= 3
y <- y[expressed,]

#Get normalized expression matrix
gexp<- y$E

gexp<- gexp[, samples]

#Get annotation
annot_chrAA<- y$genes
#write.xlsx(annot_chrAA, file = "./GSEAfiles/Annotation_for_Illumina_Probes.xlsx", rowNames=TRUE)
annot_chrAA<- cbind(PROBE_ID=rownames(annot_chrAA), annot_chrAA)
annot_chrAA<- annot_chrAA[,-3]

#Perform DEG analysis
pheno<- as.factor(c("control_sel", "control_sel", "control_sel",
                    "sel_pH647", "sel_pH647", "sel_pH647",
                    "control_acu", "control_acu", "control_acu",
                    "acu_pH68", "acu_pH68", "acu_pH68",
                    "acu_pH64", "acu_pH64", "acu_pH64"))

design <- model.matrix(~0+pheno)
colnames(design)<-gsub("pheno", "", colnames(design))

fit <- lmFit(gexp,design)
fit$genes$Symbol <- annot_chrAA$SYMBOL

#--- Contrasts of control groups against tested groups
contrasts <- makeContrasts(sel_pH647-control_sel, acu_pH68-control_acu,
                           acu_pH64-control_acu, levels=design) 

ct.fit <- eBayes(contrasts.fit(fit, contrasts))
res.fit<-decideTests(ct.fit,method="global", adjust.method="BH", p.value=0.05) #"holm", "hochberg", "bonferroni","none"

#--- Make deg table with getfitlimma
source("./Scripts/getFitLimma.R")
deg.limma_chrAA<-getFitLimma(ct.fit,res.fit)

deg.limma_chrAA$PROBEID<-rownames(deg.limma_chrAA)

#Merge deg with individual sample values from gexp
gexp<- as.data.frame(gexp)
gexp$PROBEID<-rownames(gexp)

deg.limma_chrAA2<-merge(deg.limma_chrAA, gexp, by="PROBEID")

#remove duplicates based on the probe with highest abs FC
symbols<- as.character(unique(deg.limma_chrAA2[which(duplicated(deg.limma_chrAA2$Symbol)),"Symbol"]))

result<- NULL
for (i in 1:length(symbols)) {
  gene<- symbols[i]
  mat<- deg.limma_chrAA2[which(deg.limma_chrAA2$Symbol==gene),]
  mat<-mat[order(abs(mat$`logFC.sel_pH647-control_sel`), decreasing = TRUE),]
  result<- rbind(result, mat[1,])
}

deg.limma_chrAA2<- deg.limma_chrAA2[-which(deg.limma_chrAA2$Symbol%in%symbols),]
deg.limma_chrAA2<- rbind(deg.limma_chrAA2, result)

rm(list = setdiff(ls(), c("deg.limma_chrAA2")))

####### Check how many up and down- regulated genes
length(which(deg.limma_chrAA2$`logFC.sel_pH647-control_sel`>=0.5)) # 1987
length(which(deg.limma_chrAA2$`logFC.sel_pH647-control_sel`<=(-0.5))) # 1986

#Save RDatas
save(deg.limma_chrAA2, file = "./RData/U87_AAvsNA_processedData.RData")
write.xlsx(deg.limma_chrAA2, file = "./ProcessedData/U87_AAvsNAprocessedData.xlsx")





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


