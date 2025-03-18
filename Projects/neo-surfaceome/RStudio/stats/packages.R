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
# A fast, consistent tool for working with data frame like objects, both in 
# memory and out of memory.
if (!require("TidyMultiqc")) devtools::install_github("multimeric/TidyMultiqc")
suppressPackageStartupMessages(library(TidyMultiqc))
# A package to extract and summarize data from MultiQC reports
