# Created: 2024 05 21 ; Last Modified: 2024 01 23
# KGO & MH
# ---------------------------------------------------------------------------- #
# -                         Exploratory Plots                                - #
# ---------------------------------------------------------------------------- #
setwd("E:/Lab/Collabs/Anna/Analysis2024/")


###############################################################################
#######################    Volcano plot
###############################################################################
library(openxlsx)
library(ggplot2)

# Load
Data2Plot<-read.xlsx("./ProcessedData/PANC1_AAvsNA_processedData.xlsx")

# Make a table with symbol, log2 FC and log10 Pvalue
Data2Plot<-Data2Plot[,c(4,5,6)]
colnames(Data2Plot)<-c("Symbol", "log2FoldChange", "pval")


# Define the color based on the Log2FoldChange value
Data2Plot$Color <- ifelse(Data2Plot$log2FoldChange >= 0.5, "Upregulated",
                              ifelse(Data2Plot$log2FoldChange <= (-0.5), "Downregulated", "Not-changed"))
Data2Plot$Color<-as.factor(Data2Plot$Color)

# In case you want to display specific proteins in the non-interactive plot, do:
# If multiple genes are to be highlighted
Data2Plot$Label<-ifelse(Data2Plot$Symbol %in% c("SDC1", "ADM", "HIF1A", "CA9"),
                        Data2Plot$Symbol, NA)


# Write here the title for the plot
title_plot<-"Volcano Plot PANC1 AA vs NA mRNA"

# Write here the colors for the plot
# Color names can be consulted here: https://r-graph-gallery.com/img/graph/42-colors-names.png
colors_plot<-c("darkgreen", "grey","darkred")

# Create the volcano plot (might take a while before the plot appears when opened in R due to the density of dots in the plot)
volcano<-ggplot(data=Data2Plot,aes(x=log2FoldChange,y=-log10(pval), col=Color, text=Symbol)) +
  geom_point(size = 4, alpha=0.4) + scale_color_manual(values = colors_plot)+
  ggtitle(title_plot) + geom_vline(xintercept=0.5,linetype="dashed",alpha=0.5) + xlab("log2(FC)") +
  geom_vline(xintercept=-0.5,linetype="dashed",alpha=0.5) + geom_hline(yintercept=-log10(0.05),linetype="dashed",alpha=0.5) + 
  theme_bw() + theme(axis.line=element_line(colour="black"), legend.position="top")+theme(plot.title = element_text(hjust = 0.5))+
  geom_text(aes(x=log2FoldChange,y=-log10(pval)),label=Data2Plot$Label,size=3,vjust=1.2,hjust=-0,color="black",alpha=1) 

# Save the plot (OBS: you can choose any file format, I plot as pdf whenever I want to have a vector image to open in illustrator)
ggsave(filename = "./Plots/Volcano_PANC1_AAvsNA_mRNA.png", # in here I am saving as png to have a non-vectorized image
       plot = volcano, width = 13, height = 13, units = "cm")


# Interactive Volcano - Colors in this package are not working 100&%
library(plotly)

interac_volc<-ggplotly(volcano)
interac_volc$height<-500
interac_volc$width<-700
htmltools::save_html(interac_volc, file = "./Plots/Interactive_volcano_U3034_Surface.html")

rm(list = ls())





###############################################################################
#######################    Heatmap
###############################################################################
library(gplots)
library(openxlsx)
library(heatmaply)

setwd("E:/Lab/Collabs/Anna/Analysis2024/")

# Load
Data2Plot<-read.xlsx("./ProcessedData/U87_AAvsNAprocessedData.xlsx")

# Make a table with symbol as rownames and log2 FC columns (Acute vs Chronic AA, in this case)
rownames(Data2Plot)<-Data2Plot$Symbol
Data2Plot<-Data2Plot[,c(3:5)]

# Decide which genes to display in the heatmap, makes for a better plot if the lower variance genes are removed
aux<-summary(rowMeans(abs(Data2Plot)))
threshold<- 1 # or can be aux[[5]] if you want to use genes that have higher expression than the mean
Data2Plot<- Data2Plot[-(which(rowMeans(abs(Data2Plot))<threshold)),] # Obtain the genes that pass the threshold

# Check the range of expression values to decide cap values below
range(Data2Plot)

Data2Plot[Data2Plot>=4]<-4
Data2Plot[Data2Plot<=(-4)]<-(-4)

Data2Plot<-as.matrix(Data2Plot)

pdf(file = "./Plots/Heatmap_U87AvsNA.pdf", width = 7, height = 7)
heatmap.2(Data2Plot, main = "Heatmap of U87 mRNA \n AA vs NA",
          labCol = c("Chronic", "Acute, pH 6.8", "Acute, pH 6.4"),
          labRow = FALSE, trace = "none", density.info = "none",
          dendrogram = "row", col = bluered(100), scale = "none",
          margins=c(8,4), cexCol=1)
dev.off() # run dev.off untill it appears written "null device 1" in the console

rm(list = ls())





###############################################################################
#######################    Bargraph of LogFC
###############################################################################
setwd("E:/Lab/Collabs/Anna/Analysis2024/")

library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(reshape)
library(dplyr)

# Load
Data2Plot<-read.xlsx("./ProcessedData/U87_AAvsNAprocessedData.xlsx")

# Write here the name of the genes to be plotted exactly as they appear in the expression matrix.
Genes<-c("SDC1", "ADM", "HIF1A", "CA9")

# Fish out the selected genes, and keep only the columns with symbol and lgFC information
Data2Plot<-Data2Plot[which(Data2Plot$Symbol%in%Genes), 2:5] # 2:5 is selecting Symbol and logFC columns

# Put symbols as rownames and remove the symbols column
rownames(Data2Plot)<-Data2Plot$Symbol
Data2Plot<-Data2Plot[,-1]

# Transpose the table (which is done with the t() function), to have the genes in columns and samples in rows 
Data2Plot<-as.data.frame(t(Data2Plot))

# Rename the columns to better looking names
Data2Plot$Group<-c("Chronic", "Acute_pH68", "Acute_pH64")

# Melt function rearranges the data in a way that ggplot likes
Data2Plot<-melt(Data2Plot)
colnames(Data2Plot)[2:3]<-c("Symbol", "Log2FC")

# Plotting and saving
p<-ggplot(Data2Plot, aes(x=Group, y=Log2FC, fill=Group)) +
  geom_bar(stat = "identity", width=0.6, color="black", linewidth = 0.4) +
  theme_bw()+ scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  ggtitle("Genes of interest")+ scale_fill_brewer(palette="Blues")+
  scale_x_discrete(guide = guide_axis(angle = 45)) + facet_wrap(~Symbol, scales = "free")

ggsave("./Plots/Specific_Genes.png", plot = p, width = 15, height = 15, units = "cm")

rm(list = ls())





###############################################################################
##########    Bargraph of individual abundances per sample in linear scale
###############################################################################
# Load
Data2Plot<-read.xlsx("./ProcessedData/U87_AAvsNAprocessedData.xlsx")

# Write here the name of the genes to be plotted exactly as they appear in the expression matrix.
Genes<-c("SDC1", "ADM", "HIF1A", "CA9")

# Fish out the selected genes, and keep the columns with symbol and sample-wise expression
Data2Plot<-Data2Plot[which(Data2Plot$Symbol%in%Genes), c(2,15:29)]

# Put symbols as rownames and remove the symbols column
rownames(Data2Plot)<-Data2Plot$Symbol
Data2Plot<-Data2Plot[,-1]

# Transpose the table (which is done with the t() function), to have the genes in columns and samples in rows 
Data2Plot<-as.data.frame(t(Data2Plot))

# Provide the name of the group each sample belongs to
Data2Plot$Group<- factor(c("control_sel", "control_sel", "control_sel",
                              "sel_pH647", "sel_pH647", "sel_pH647",
                              "control_acu", "control_acu", "control_acu",
                              "acu_pH68", "acu_pH68", "acu_pH68",
                              "acu_pH64", "acu_pH64", "acu_pH64"),
                            levels=c("control_acu", "acu_pH64", "acu_pH68",
                                     "control_sel", "sel_pH647"))

# Transform the normalized log2 values into linear scale
Data2Plot[,1:4]<-2^Data2Plot[,1:4]

# Melt function rearranges the data in a way that ggplot likes
Data2Plot<-melt(Data2Plot)
colnames(Data2Plot)[2:3]<-c("Symbol", "Lin")

# This part obtains the mean and standard deviation for each gene
df<-Data2Plot %>% group_by(Group, Symbol) %>% summarize(mean(Lin),sd(Lin))
colnames(df)[3:4]<-c("mean_lin", "sd_lin")

# Colors for the plot are defined below with the scale_fill_brewer function.
# Find color paletes here: https://rstudio-pubs-static.s3.amazonaws.com/5312_98fc1aba2d5740dd849a5ab797cc2c8d.html

# Plot and save
p<-ggplot(df, aes(x=Group, y=mean_lin, fill=Group)) +
  geom_bar(stat = "identity", width=0.6, color="black", linewidth = 0.4) +
  theme_bw()+ scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  ggtitle("Genes of interest")+ scale_fill_brewer(palette="Blues")+
  scale_x_discrete(guide = guide_axis(angle = 45)) + facet_wrap(~Symbol, scales = "free")+
  geom_errorbar(aes(ymin=mean_lin-sd_lin, ymax=mean_lin+sd_lin), width=.2,
                position=position_dodge(.3))

ggsave("./Plots/Specific_Genes_Linearscale_individualSamples.png", plot = p, width = 15, height = 15, units = "cm")

rm(list = ls())





###############################################################################
#######################    Gene signature score
###############################################################################
library(hacksig)

# Load data (OBS: for IvyGap, load the RData)
load("./ProcessedData/IvyGap_processedData.RData")
rm(IvyGap, clin_IvyGap)
rownames(sel_IvyGap)<-sel_IvyGap$Sample

# Define which genes are part of the gene signature
Lipid_droplet_signature<- c("HILPDA", "PLIN2", "PLIN3", "G0S2", "PLIN1", "DGAT1", "SOAT1")


# Use transpose function to put samples in columns and genes in rows
IvyGApforSig<-as.data.frame(t(sel_IvyGap[,18:ncol(sel_IvyGap)]))


# Perform the score calculation
Lipid_score<-hack_sig(IvyGApforSig, signatures = list(Lipid_droplet_signature), method = "zscore")
colnames(Lipid_score)<-c("Sample", "LDsig_score")

# Add histology information next to the LD signature score
Lipid_score<-merge(Prim_IG[,1:2], Lipid_score, by="Sample")

# Plot and save
p <- ggplot(Lipid_score, aes(x=Histology, y=LDsig_score)) + 
  geom_boxplot(na.rm = TRUE, size=1, width=0.5, color="black", linewidth = 0.4)+
  labs(title="Lipid droplet signature scores in IvyGap data \n (pval = T.test, ref group = PPC)")+
  theme_bw()+ scale_x_discrete(guide = guide_axis(angle = 45))+
  stat_compare_means(label = "p.signif",
                     method = "t.test", ref.group = "Pseudopalisading cells") # Add pairwise comparisons p-value

ggsave(p, filename = "./Plots/LipidDroplet_SignatureScore_IvyGap.png",
       width = 13, height = 15, units = "cm", dpi = 600)

rm(list = ls())
