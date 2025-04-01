if (!requireNamespace("RSQLite", quietly = TRUE))
   BiocManager::install("RSQLite", update = FALSE)
suppressPackageStartupMessages(library(RSQLite))

if (!requireNamespace("DBI", quietly = TRUE))
  BiocManager::install("DBI", update = FALSE)
suppressPackageStartupMessages(library(DBI))

if (!requireNamespace("AnnotationDbi", quietly = TRUE))
  BiocManager::install("AnnotationDbi", update = FALSE)
suppressPackageStartupMessages(library(AnnotationDbi))

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db", update = FALSE)
suppressPackageStartupMessages(library(org.Hs.eg.db))

if (!requireNamespace("rbioapi", quietly = TRUE))
  install.packages("rbioapi")
suppressPackageStartupMessages(library(rbioapi))

if (!requireNamespace("PFAM.db", quietly = TRUE))
  BiocManager::install("PFAM.db", update = FALSE)
suppressPackageStartupMessages(library(PFAM.db))

if (!requireNamespace("tidyverse", quietly = TRUE))
  BiocManager::install("tidyverse", update = FALSE)
suppressPackageStartupMessages(library(tidyverse))