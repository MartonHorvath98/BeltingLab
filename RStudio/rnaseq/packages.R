################################################################################
# Load required packages                                                       #
################################################################################

# 1.) Loader packages ##########################################################
if (!require("BiocManager")) install.packages("BiocManager")
suppressPackageStartupMessages(library(BiocManager))
# Provides tools for managing Bioconductor packages
if (!require("openxlsx")) install.packages("openxlsx")
suppressPackageStartupMessages(library(openxlsx)) 
# Simplifies the creation, writing, and editing of '.xlsx' files
if (!require("devtools")) install.packages("devtools")
suppressPackageStartupMessages(library(devtools)) 
# Collection of package development tools.

# 2.) Utility packages #########################################################

if (!require("crayon")) install.packages("crayon")
suppressPackageStartupMessages(library(crayon))
# Provides a set of functions to style console output
if (!require("dplyr")) install.packages("dplyr")
suppressPackageStartupMessages(library(dplyr))
# tool for working with data frame like objects
if (!require("tidyr")) install.packages("tidyr")
suppressPackageStartupMessages(library(tidyr))
# contains tools for changing the shape (pivoting) and hierarchy (nesting and 
# 'unnesting') of a dataset
if (!require("tibble")) install.packages("tibble")
suppressPackageStartupMessages(library(tibble))
# provides a 'tbl_df' class that offers better formatting than the default data 
# frame
if (!require("plyr")) install.packages("plyr")
suppressPackageStartupMessages(library(plyr))
# provides a set of tools for splitting, applying, and combining data
if (!require("stringr")) install.packages("stringr")
suppressPackageStartupMessages(library(stringr))
# provides a consistent, simple and easy to use set of wrappers around the
# 'stringi' package
if (!require("forcats")) install.packages("forcats")
suppressPackageStartupMessages(library(forcats))
# provides a set of tools for working with factors
if (!require("numbers")) install.packages("numbers")
suppressPackageStartupMessages(library(numbers))
# provides a set of functions to format numbers
if (!require("conicfit")) install.packages("conicfit")
suppressPackageStartupMessages(library(conicfit))
# provides a set of functions to fit a conic section to a set of points

# 3.) Differential expression ##################################################

if (!require("DESeq2")) BiocManager::install("DESeq2", update = F)
suppressPackageStartupMessages(library(DESeq2))
# provides methods to test for differential expression in count data
if (!require("fdrtool")) install.packages("fdrtool")
suppressPackageStartupMessages(library(fdrtool))
# provides a set of functions to estimate the false discovery rate (FDR) and
# q-values from a list of p-values
if (!require("edgeR")) install.packages("edgeR")
suppressPackageStartupMessages(library(edgeR))
# provides methods to test for differential expression in count data
if (!require("preprocessCore")) BiocManager::install("preprocessCore", update = F)
suppressPackageStartupMessages(library(preprocessCore))
# Provides a set of functions to preprocess microarray data
if (!require("WGCNA")) install.packages("WGCNA")
suppressPackageStartupMessages(library(WGCNA))
# provides a set of functions to perform weighted gene co-expression network
# analysis
if (!require("GenomicFeatures")) BiocManager::install("GenomicFeatures", update = F)
suppressPackageStartupMessages(library(GenomicFeatures))
# provides a set of functions to work with genomic features
if (!require("granulator")) BiocManager::install("granulator", update=F)
suppressPackageStartupMessages(library(granulator))
# provides a set of functions to perform gene set enrichment analysis

# 4.) Gene ontology enrichment ################################################
if (!require("clusterProfiler")) install.packages("clusterProfiler")
suppressPackageStartupMessages(library(clusterProfiler))
# provides a set of functions to perform gene set enrichment analysis
if (!require("org.Hs.eg.db")) install.packages("org.Hs.eg.db")
suppressPackageStartupMessages(library(org.Hs.eg.db))
# provides a set of functions to map between gene identifiers and gene symbols
if (!require("GO.db")) install.packages("GO.db")
suppressPackageStartupMessages(library(GO.db))
# provides a set of functions to map between gene identifiers and gene ontology
# terms
if (!require("topGO")) BiocManager::install("topGO", update = F)
suppressPackageStartupMessages(library(topGO))
# provides a set of functions to perform gene ontology enrichment analysis
if (!require("rrvgo")) BiocManager::install("rrvgo", update = F)
suppressPackageStartupMessages(library(rrvgo))
# provides a set of functions to perform gene ontology enrichment analysis
if (!require("GOSemSim")) BiocManager::install("GOSemSim", update = F)
suppressPackageStartupMessages(library(GOSemSim))
# provides a set of functions to calculate semantic similarity between gene
# ontology terms

# 5.) Data visualization ######################################################

if (!require("ggplot2")) install.packages("ggplot2")
suppressPackageStartupMessages(library(ggplot2))
# provides a set of functions to create plots
if (!require("ggbiplot")) install.packages("ggbiplot")
suppressPackageStartupMessages(library(ggbiplot))
# provides a set of functions to create biplots
if (!require("pheatmap")) install.packages("pheatmap")
suppressPackageStartupMessages(library(pheatmap))
# provides a set of functions to create heatmaps
if (!require("RColorBrewer")) install.packages("RColorBrewer")
suppressPackageStartupMessages(library(RColorBrewer))
# provides a set of functions to create color palettes
if (!require("ggpubr")) install.packages("ggpubr")
suppressPackageStartupMessages(library(ggpubr))
# provides a set of functions to create publication-ready plots
if (!require("ggrepel")) install.packages("ggrepel")
suppressPackageStartupMessages(library(ggrepel))
# provides a set of functions to create plots with text labels
if (!require("GOplot")) install.packages("GOplot")
suppressPackageStartupMessages(library(GOplot))
# provides a set of functions to create gene ontology plots
if (!require("VennDiagram")) install.packages("VennDiagram")
suppressPackageStartupMessages(library(VennDiagram))
# provides a set of functions to create Venn diagrams
if (!require("ComplexUpset")) install.packages("ComplexUpset", 
                                               repos = BiocManager::repositories(),
                                               type = "source")
suppressPackageStartupMessages(library(ComplexUpset))
# provides a set of functions to create complex UpSet plots
if (!require("ComplexHeatmap")) install.packages("ComplexHeatmap", 
                                                 repos = BiocManager::repositories(),
                                                 type = "source")
suppressPackageStartupMessages(library(ComplexHeatmap))
# provides a set of functions to create complex heat maps
if (!require("circlize")) install.packages("circlize")
suppressPackageStartupMessages(library(circlize))
# provides a set of functions to create circular plots
if (!require("TidyMultiqc")) install.packages("TidyMultiqc", 
                                               repos = BiocManager::repositories(),
                                               type = "source")
suppressPackageStartupMessages(library(TidyMultiqc))
# provides a set of functions to create MultiQC plots
if (!require("igraph")) install.packages("igraph")
suppressPackageStartupMessages(library(igraph))
# provides a set of functions to create network plots
