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
if (!require("stringr")) install.packages("stringr")
suppressPackageStartupMessages(library(stringr))
# provides a consistent, simple and easy to use set of wrappers around the
# 'stringi' package 
if (!require("DESeq2")) install.packages("DESeq2")
suppressPackageStartupMessages(library(DESeq2))
# provides methods to test for differential expression in count data
if (!require("fdrtool")) install.packages("fdrtool")
suppressPackageStartupMessages(library(fdrtool))
# provides a set of functions to estimate the false discovery rate (FDR) and
# q-values from a list of p-values
if (!require("edgeR")) install.packages("edgeR")
suppressPackageStartupMessages(library(edgeR))
# provides methods to test for differential expression in count data
if (!require("ggplot2")) install.packages("ggplot2")
suppressPackageStartupMessages(library(ggplot2))
# provides a set of functions to create plots
if (!require("ggbiplot")) install.packages("ggbiplot")
suppressPackageStartupMessages(library(ggbiplot))
# provides a set of functions to create biplots
