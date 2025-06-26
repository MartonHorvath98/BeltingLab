# 1.  Transcriptomic Analysis of the effects of acidosis on Glioblastoma cell lines and patient samples

- [1.  Transcriptomic Analysis of the effects of acidosis on Glioblastoma cell lines and patient samples](#1--transcriptomic-analysis-of-the-effects-of-acidosis-on-glioblastoma-cell-lines-and-patient-samples)
  - [1.1. Input data](#11-input-data)
  - [1.2. Directory tree](#12-directory-tree)
- [2. Bioinformatical pipeline](#2-bioinformatical-pipeline)
  - [2.1. Data preprocessing \& differential gene expression (DGE) analysis](#21-data-preprocessing--differential-gene-expression-dge-analysis)
      - [Overview:](#overview)
      - [Main steps:](#main-steps)
      - [Input data:](#input-data)
      - [Output:](#output)
  - [2.2. DGE visualization and pathway analysis](#22-dge-visualization-and-pathway-analysis)
      - [Overview:](#overview-1)
      - [Main Steps:](#main-steps-1)
      - [Input data:](#input-data-1)
      - [Output:](#output-1)
  - [2.3. Addition GSEA visualization](#23-addition-gsea-visualization)
      - [Overview:](#overview-2)
      - [Main Steps:](#main-steps-2)
      - [Input data:](#input-data-2)
      - [Output:](#output-2)
  - [2.4. GSEA clustering and network analysis](#24-gsea-clustering-and-network-analysis)
      - [Overview:](#overview-3)
      - [Main Steps:](#main-steps-3)
      - [Input data:](#input-data-3)
      - [Output:](#output-3)
  - [2.5. Gene signature analysis](#25-gene-signature-analysis)
      - [Overview:](#overview-4)
      - [Main Steps:](#main-steps-4)
      - [Inputs:](#inputs)
      - [Outputs:](#outputs)

## 1.1. Input data

Three main, different data sources were used in these analyses:
1. **Affymetrix Clariom D Pico Gene Array** - performed on pooled lasermicrodissected samples from lipid droplet rich (n=5) and lipid droplet lacking (n=5) matched pathological samples (hereinafter called **CCLD**); three patient-derived primary cell lines grown as 3D vs. 2D cultures (**HGCC**); and a pancreatic cancer cell line, PANC1, stimulated by a low-pH environment (pH6.4) for 10 weeks (adapted acidosis, AA) (**PANC1**).
2. **Illumina HumanHT-12 v4 Expression BeadChip** - used for gene expression analysis of U87 cells grown either in low-pH or hypoxic growth conditions (**U87**).

> [!NOTE]
> *All raw and processed data used in this analysis have been uploaded to NCBI's public data repository, the Gene Expression Omnibus (GEO).*

The data are available under the unique IDs **GSE300758** (**CCLD**), **GSE300765** (**U87**), **GSE300768** (**PANC1**), and **GSE300771** (**HGCC**). The accession of the datasets is scheduled to stay private until 30th September, 2025, unless the manuscript enters review earlier.

## 1.2. Directory tree

```bash
.
├── README.md
├── RStudio/
│   ├── RData/
│   ├── Results/
│   ├── DEGanalysis.R
│   ├── GSEAvenn.R
│   ├── SignatureAnalysis.R
│   ├── clusterTerms.R
│   ├── functions.R
│   ├── packages.R
│   └── preprocessing.R
└── data/
    ├── processed/
    └── raw/
```

# 2. Bioinformatical pipeline 

## 2.1. Data preprocessing & differential gene expression (DGE) analysis

#### Overview:
The script [`preprocessing.R`](RStudio/preprocessing.R) performs preprocessing, normalization, and differential expression analysis (DEA) for multiple transcriptomic datasets from Affymetrix Clariom D Human Pico and Illumina BeadChip platforms.

#### Main steps:
1. Environment Setup:
     - Sets working directory and loads required packages and custom functions.

2. DGE Analysis each dataset:
     - Svenja's CC +LD and CC no LD (Affymetrix)
     - Hugo's Primary cells 2D vs 3D (Affymetrix)
     - U87 Chronic Acidosis AA vs NA & HOX vs NOX (Illumina)
     - PANC1 Chronic Acidosis AA vs NA (Affymetrix)

For each dataset:
  - Creates distinct results (tables) directories.
  - Loads raw files.
  - Normalize raw data and perform transcript filtering.
  - Conduct differential expression analysis using limma.
  - Remove duplicate gene symbols/IDs based on log2 fold change or t-stat.
  - Summarize up- and down-regulated genes for each comparison.
  - Save processed data as RData and Excel files for downstream analysis.

Requirements:
  - R packages: oligo, affycoretools, limma, stats
  - Annotations: clariomdhumantranscriptcluster.db (v8.8.0), illuminaHumanv4.db (v1.26.0), org.Hs.eg.db (v3.20.0)
  - Custom functions (from [`functions.R`](RStudio/functions.R)): normalizeTranscript, normalizeIllumina, limmaDEA, removeDuplicates

#### Input data:
  - Raw `.cel` files for Affirmetrix experiments, and calculated probe and control intensity summaries for the Illumina experiment

#### Output:
  - Processed RData and Excel files for each dataset and comparison

## 2.2. DGE visualization and pathway analysis

#### Overview:
The script [`DEGanalysis.R`](RStudio/DEGanalysis.R) performs the visual analysis of the DGE analysis results and downstream pathway enrichment for multiple transcriptomics datasets (Affymetrix Clariom D, Illumina BeadChip). This script performs a comprehensive analysis including differential gene expression (DGE), over-representation analysis (ORA), gene set enrichment analysis (GSEA), and visualization of results. The workflow is modular and supports multiple datasets and experimental designs.

#### Main Steps:
1. Environment Setup:
     - Sets working directory and loads required packages and custom functions.
     - Loads MSigDB gene set databases (Hallmark, GO:BP, KEGG, Reactome).

2. Pathway Analysis for each dataset:
     - Svenja's CC +LD and CC no LD (Affymetrix)
     - Hugo's Primary cells 2D vs 3D (Affymetrix)
     - U87 Chronic Acidosis AA vs NA & HOX vs NOX (Illumina)
     - PANC1 Chronic Acidosis AA vs NA (Affymetrix)

For each dataset:
  - Creates distinct results (tables and plots) directories
  - Prepares metadata and sample annotation.
  - Identifies and annotates DEGs, and perform exploratory visual analysis (PCA, Venn, Volcano plots).
  - Extracts gene lists for enrichment analysis and runs Over-Representation Analysis (ORA) and Gene Set Enrichment Analysis (GSEA) using MSigDB gene sets.

Requirements:
  - R packages: dplyr, tidyr, stringr, openxlsx, msigdbr, org.Hs.eg.db, cowplot, ggplot2, etc.
  - Custom functions (from [`functions.R`](RStudio/functions.R)): get_significance, get_genelist, run_ora, extract_ora_results, run_gsea, extract_gsea_results, plot_pca, plot_venn, plot_vulcan, getEnrichmentTable, plotRunningScore, plotGeneRank, plotClusters.

#### Input data:
  - Expression matrices and DEG tables for each dataset.

#### Output:
  - Differential expression tables, pathway enrichment results, and publication-ready plots.

## 2.3. Addition GSEA visualization

#### Overview:
Besides the [`DEGanalysis.R`](RStudio/DEGanalysis.R) script, [`GSEAvenn.R`](RStudio/GSEAvenn.R) script performs additional visual analysis of the pathway enrichment results, creating Venn diagrams and heatmap across multiple experimental conditions. The Venn diagrams and heatmaps visualize the overlap and trends of enriched GO terms and pathways, as well as the expression of selected genes of interest.

#### Main Steps:
1. Comparison of Uppsala cells and Lipid loaded cells:
     - Merges GSEA and GO enrichment results for four conditions (U3017, U3047, U3054, CCLD).
     - Identifies overlapping and unique gene sets across conditions.
     - Assigns each gene set to a Venn diagram region based on its presence in the conditions.
     - Visualizes the overlap using a custom Venn diagram with ggplot2, ggforce, and scatterpie.
     - Annotates each region with the number of gene sets and their regulation trends (UP, DOWN, CHANGE).
     - Saves the Venn diagram and exports tables of gene sets for each region.
     - Generates a heatmap of selected genes of interest, grouped by functional category.

2. Comparison of AA vs CA vs Hypoxia in U87 cells:
    - Merges GSEA and GO enrichment results for four U87 cell conditions (pH6.4_CA, pH6.4_AA, pH6.8_AA, Hypoxia).
    - Identifies overlapping and unique gene sets across conditions.
    - Assigns each gene set to a Venn diagram region based on its presence in the conditions.
    - Visualizes the overlap using a custom Venn diagram with ggplot2, ggforce, and scatterpie.
    - Annotates each region with the number of gene sets and their regulation trends (UP, DOWN, CHANGE).
    - Saves the Venn diagram and exports tables of gene sets for each region.

Requirements:
  - R packages: dplyr, tidyr, tibble, ggplot2, ggforce, scatterpie, ComplexHeatmap etc.

#### Input data:
  - GSEA and GO enrichment result objects for each condition.
  - List of genes of interest with functional categories.

#### Output:
  - Venn diagrams (PNG, SVG) showing overlap of enriched gene sets.
  - Heatmap of selected genes of interest.
  - Excel tables of gene sets for each Venn region.

## 2.4. GSEA clustering and network analysis

#### Overview:

The [`clusterTerms.R`](RStudio/clusterTerms.R) script performs clustering and visualization of enriched pathways,  terms, and hallmarks based on gene set enrichment analysis (GSEA) results for U87 and PANC1 cell lines under various conditions (acidosis, hypoxia).

#### Main Steps:

1. Calculate Cohen's similarity matrix between all gene sets (GO terms, MSigDb hallmarks, KEGG, Reactome pathways).
2. Cluster terms and pathways based on the similarity matrix for each experimental comparison.  
3. Identify representative terms for each cluster.
4. Visualize clusters as networks and circos plots, highlighting selected pathways and genes of interest.  
5. Generate volcano plots for selected genes.  

Requirements:
  - R packages: igraph, RCy3, etc.
  - custom functions: cohen_kappa, get_cluster, get_cluster_representative, filter_graph, plot_network, getCircplotData, plotCircplot, plot_vulcan

#### Input data:
   - GSEA and GO enrichment results for U87 and PANC1 cell lines
   - Differential expression results for each corresponding comparison 
   - Predefined lists of pathways/terms and genes of interest 

#### Output:
  - Clustered enrichment tables (xlsx)
  - Network visualizations (PNG, SVG)
  - Volcano plots for genes of interest (PNG) 

## 2.5. Gene signature analysis

#### Overview:
The [`SignatureAnalysis.R`] script performs gene signature analysis for lipid droplet genes across multiple glioblastoma datasets (IvyGap, TCGA, CGGA). It calculates gene signature scores using `zscore`, assigns 
molecular subtypes using `ssgsea`, and generates visualizations and statistical summaries.

The cohorts used for this segment of the project include RNA-seq data from all three datasets. In the case of Ivygap and CGGA the normalized count reads from the pre-processed data (sequence alignment and transcript 
abundance estimation) were log2 transformed after adding a 0.5 pseudocount (to avoid infinite value upon log transformation). 

In the case of the IvyGap dataset primary tumors were filtered based on 'Histology status' and samples that have location stated (Ivygap does not have IDH status information). This way out of the 270 samples, 122 were retained as primary GBM. In case of the TCGA dataset outof the 548 samples 142 primary GBM was kept, and in the CGGA dataset out of the 1018 samples 183.

#### Main Steps:

1. Load the gene signature and the GBM transcriptomic classification (Verhaat et al., hallmark gene sets).
2. Load, preprocess, and filter IvyGap, TCGA, and CGGA datasets for primary GBM samples.
3. Calculate lipid droplet gene signature scores using z-score method.
4. Assign molecular subtypes using ssGSEA and best-fit approaches.
5. Generate and save plots for signature scores by histology and subtype.
6. Export statistical test results and processed data to Excel files.

Requirements:
  - R packages: hacksig, ggplot, msigdbr etc.
  - commercial functions: hack_sig, and custom functions: score_plot

#### Inputs:
   - data/gene-signature-new.txt: List of lipid droplet genes.
   - msigdbr_df: Data frame of hallmark gene sets.
   - IvyGap, TCGA, CGGA expression and clinical data files.

#### Outputs:
   - PNG and SVG plots of signature scores by histology and subtype.
   - Excel files with statistical test results and processed datasets.
   - RData files with processed data objects.