# Environment: R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("BiocParallel", quietly = TRUE)) BiocManager::install("BiocParallel", update = F)
suppressPackageStartupMessages(library(BiocParallel))
# Base utility packages
if (!requireNamespace("dplyr", quietly = TRUE)) BiocManager::install("dplyr", update = F)
suppressPackageStartupMessages(library(dplyr))
if (!requireNamespace("tidyr", quietly = TRUE)) BiocManager::install("tidyr", update = F)
suppressPackageStartupMessages(library(tidyr))
if (!requireNamespace("stringr", quietly = TRUE)) BiocManager::install("stringr", update = F)
suppressPackageStartupMessages(library(stringr))
if (!requireNamespace("openxlsx", quietly = TRUE)) BiocManager::install("openxlsx", update = F)
suppressPackageStartupMessages(library(openxlsx))
#VCF manipulation
if (!requireNamespace("vcfR", quietly = TRUE)) BiocManager::install("vcfR", update = F)
suppressPackageStartupMessages(library(vcfR))
if (!requireNamespace("maftools", quietly = TRUE)) BiocManager::install("maftools", update = F)
suppressPackageStartupMessages(library(maftools))
#Annotation
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi", update = F)
suppressPackageStartupMessages(library(AnnotationDbi))
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db", update = F)
suppressPackageStartupMessages(library(org.Hs.eg.db))
#Isoform analysis
if (!requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE)) BiocManager::install("IsoformSwitchAnalyzeR", update = F)
suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))
if (!requireNamespace("ballgown", quietly=TRUE)) BiocManager::install("ballgown", update = F)
suppressPackageStartupMessages(library(ballgown))
#Visualization
if (!requireNamespace("ggplot2", quietly = TRUE)) BiocManager::install("ggplot2", update = F)
suppressPackageStartupMessages(library(ggplot2))
if (!requireNamespace("cowplot", quietly = TRUE)) BiocManager::install("cowplot", update = F)
suppressPackageStartupMessages(library(cowplot))
if (!requireNamespace("ggpubr", quietly = TRUE)) BiocManager::install("ggpubr", update = F)
suppressPackageStartupMessages(library(ggpubr))

