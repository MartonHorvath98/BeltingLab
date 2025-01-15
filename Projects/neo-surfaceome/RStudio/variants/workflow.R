# ---------------------------------------------------------------------------- #
# -            SET UP WORKING DIRECTORY AND DIRECTORY TREE                   - #
# ---------------------------------------------------------------------------- #
# 1.) Set downstream path
## Create the sub folders for: results, data, and pictures
data <- "data"
if (!dir.exists(file.path(data))) {
  dir.create(file.path(data)) # create the data folder
} 
plots_dir <- "plots" # plots directory
results_dir <- "results" # results directory
date <- format(Sys.Date(), "%Y-%m-%d") # get the current date
if (!dir.exists(file.path(date))) {
  dir.create(file.path(date)) # create the dated results folder
  dir.create(file.path(date, results_dir)) # create the results folder
  dir.create(file.path(date, plots_dir)) # create the plots folder
}

# 2.) Load functions for the analyses and required packages
source("functions.R")
source("packages.R")
source("helpers.R")

# 3.) Load analysis data              
# ------- 1) generated with BCFtools -----------------
bcftools.path <- file.path("../../results/05_variant/BCFtools/")
bcftools.files <- list.files(bcftools.path, pattern = ".vcf.gz$", full.names = T)
# exclude hypoxia samples
bcftools.files <- bcftools.files[!grepl("2DH", bcftools.files)]
bcftools.data <- lapply(bcftools.files, function(x) { 
  raw <- read.vcfR(x)
  tidy <- vcfR2tidy(raw)
  return(tidy)
})
names(bcftools.data) <- conditions
# Format the data
bcftools.table <- merge_vcf(bcftools.data)
summary(bcftools.table)
 
# ------- 2) generated with GATK -----------------
gatk.path <- file.path("../../results/05_variant/GATK/")
gatk.files <- list.files(gatk.path, pattern = ".vcf.gz$", full.names = T)
# exclude hypoxia samples
gatk.files <- gatk.files[!grepl("2DH", gatk.files)]
gatk.data <- lapply(gatk.files, function(x) { 
  raw <- read.vcfR(x)
  tidy <- vcfR2tidy(raw)
  return(tidy)
})
names(gatk.data) <- conditions
# Format the data
gatk.table <- merge_vcf(gatk.data)
summary(gatk.table)

# 4.) Create a complete data table and visualize basic metrics
vcf.table <- rbind(bcftools.table %>% dplyr::mutate(method = "bcftools"),
                    gatk.table %>% dplyr::mutate(method = "gatk")) %>% 
  dplyr::mutate(method = factor(method, levels = c("bcftools","gatk"))) %>% 
  dplyr::filter(type == "SNP")

vcf.summary <- summarise_vcf(vcf.table, conv, conv.class)

vcf.summary$plot <- plot_summary(vcf.summary$res)

# Save the plot
ggsave(plot = vcf.summary$plot, filename = file.path(date, plots_dir, "variant_summary.png"),
       device = "png", width = 14, height = 8, units = "in", dpi = 300)

# 5.) Filter the data
### Load the SURFME data
SURFME <- read.xlsx(file.path("data","SURFME_v2023.xlsx"))
colnames(SURFME) <- c("UniProtID", "Description","Emtry", "geneSymbol",
                      "EnsemblID", "GPI", "SinglePass", "MultiPass", "Cellulare_membrane",
                      "Extracellular_domain", "GOCellSurf", "GOExternal", "Surfy",
                      "GOPlasmaMembrane")
SURFME_filter <-  unique(SURFME$Ensembl, na.rm = T)
### Load the somatic data
VI_3429_593_somatic <- load_wgs(file.path("data","593-tumor-tissue-consensus.tsv"))
VI_3429_673_somatic <- load_wgs(file.path("data","673-tumor-tissue-consensus.tsv"))

vcf.filtered.593 <- vcf.table %>% 
  dplyr::filter(sample %in% conditions[grepl("593", conditions)]) %>%
  dplyr::group_split(sample) %>% 
  setNames(., conditions[grepl("593", conditions)])

vcf.filtered.673 <- vcf.table %>%
  dplyr::filter(sample %in% conditions[grepl("673", conditions)]) %>%
  dplyr::group_split(sample) %>% 
  setNames(., conditions[grepl("673", conditions)])


vcf.filtered.593 <- lapply(vcf.filtered.593, filter_vcf, VI_3429_593_somatic, SURFME_filter)
vcf.filtered.673 <- lapply(vcf.filtered.673, filter_vcf, VI_3429_673_somatic, SURFME_filter)

vcf.filter.summary.593 <- lapply(vcf.filtered.593, function(x){
  lapply(x, function(y){
    return(summarise_vcf(y, conv, conv.class)$res)})
})


# 4.) Calculate the number of variants after every filtering step
variant_num <- data.frame(
  "total" = c(sapply(vcf.filtered.593, function(x) {x[["df"]] %>% dplyr::pull(variant_ID) %>% unique() %>% length()}),
              sapply(vcf.filtered.673, function(x) {x[["df"]] %>% dplyr::pull(variant_ID) %>% unique() %>% length()})),
  "total.bcftools" = c(sapply(vcf.filtered.593, function(x) {x[["df"]] %>% 
                           dplyr::filter(method == "bcftools") %>% 
                           dplyr::pull(variant_ID) %>% unique() %>% length()}),
                       sapply(vcf.filtered.673, function(x) {x[["df"]] %>% 
                           dplyr::filter(method == "bcftools") %>% 
                           dplyr::pull(variant_ID) %>% unique() %>% length()})),
  "total.gatk" = c(sapply(vcf.filtered.593, function(x) {x[["df"]]  %>% 
                           dplyr::filter(method == "gatk") %>% 
                           dplyr::pull(variant_ID) %>% unique() %>% length()}),
                      sapply(vcf.filtered.673, function(x) {x[["df"]] %>% 
                          dplyr::filter(method == "gatk") %>% 
                          dplyr::pull(variant_ID) %>% unique() %>% length()})),
  "total.shared" = c(sapply(vcf.filtered.593, function(x) {x[["df"]] %>% 
                         dplyr::group_by(variant_ID) %>% 
                         dplyr::filter(n() > 1) %>% 
                         dplyr::pull(variant_ID) %>% unique() %>% length()}),
                     sapply(vcf.filtered.673, function(x) {x[["df"]] %>% 
                         dplyr::group_by(variant_ID) %>% 
                         dplyr::filter(n() > 1) %>% 
                         dplyr::pull(variant_ID) %>% unique() %>% length()})),
  "somatic" = c(sapply(vcf.filtered.593, function(x) {x[["somatic"]] %>% dplyr::pull(variant_ID) %>% unique() %>% length()}),
              sapply(vcf.filtered.673, function(x) {x[["somatic"]] %>% dplyr::pull(variant_ID) %>% unique() %>% length()})),
  "somatic.bcftools" = c(sapply(vcf.filtered.593, function(x) {x[["somatic"]] %>% 
                        dplyr::filter(method == "bcftools") %>% 
                        dplyr::pull(variant_ID) %>% unique() %>% length()}),
                    sapply(vcf.filtered.673, function(x) {x[["somatic"]] %>% 
                        dplyr::filter(method == "bcftools") %>% 
                        dplyr::pull(variant_ID) %>% unique() %>% length()})),
  "somatic.gatk" = c(sapply(vcf.filtered.593, function(x) {x[["somatic"]]  %>% 
                        dplyr::filter(method == "gatk") %>% 
                        dplyr::pull(variant_ID) %>% unique() %>% length()}),
                    sapply(vcf.filtered.673, function(x) {x[["somatic"]] %>% 
                        dplyr::filter(method == "gatk") %>% 
                        dplyr::pull(variant_ID) %>% unique() %>% length()})),
  "somatic.shared" = c(sapply(vcf.filtered.593, function(x) {x[["somatic"]] %>% 
                        dplyr::group_by(variant_ID) %>% 
                        dplyr::filter(n() > 1) %>% 
                        dplyr::pull(variant_ID) %>% unique() %>% length()}),
                    sapply(vcf.filtered.673, function(x) {x[["somatic"]] %>% 
                        dplyr::group_by(variant_ID) %>% 
                        dplyr::filter(n() > 1) %>% 
                        dplyr::pull(variant_ID) %>% unique() %>% length()})),
  "surface" = c(sapply(vcf.filtered.593, function(x) {x[["surface"]] %>% dplyr::pull(variant_ID) %>% unique() %>% length()}),
              sapply(vcf.filtered.673, function(x) {x[["surface"]] %>% dplyr::pull(variant_ID) %>% unique() %>% length()})),
  "surface.bcftools" = c(sapply(vcf.filtered.593, function(x) {x[["surface"]] %>% 
      dplyr::filter(method == "bcftools") %>% 
      dplyr::pull(variant_ID) %>% unique() %>% length()}),
      sapply(vcf.filtered.673, function(x) {x[["surface"]] %>% 
          dplyr::filter(method == "bcftools") %>% 
          dplyr::pull(variant_ID) %>% unique() %>% length()})),
  "surface.gatk" = c(sapply(vcf.filtered.593, function(x) {x[["surface"]]  %>% 
      dplyr::filter(method == "gatk") %>% 
      dplyr::pull(variant_ID) %>% unique() %>% length()}),
      sapply(vcf.filtered.673, function(x) {x[["surface"]] %>% 
          dplyr::filter(method == "gatk") %>% 
          dplyr::pull(variant_ID) %>% unique() %>% length()})),
  "surface.shared" = c(sapply(vcf.filtered.593, function(x) {x[["surface"]] %>% 
      dplyr::group_by(variant_ID) %>% 
      dplyr::filter(n() > 1) %>% 
      dplyr::pull(variant_ID) %>% unique() %>% length()}),
      sapply(vcf.filtered.673, function(x) {x[["surface"]] %>% 
          dplyr::group_by(variant_ID) %>% 
          dplyr::filter(n() > 1) %>% 
          dplyr::pull(variant_ID) %>% unique() %>% length()})))
      


# Plot the top 10 mutated genes in both patients
(top_genes.593 <- plot_top_genes(vcf.filtered.593, VI_3429_593_somatic, n = 10, "Patient 1"))
top_genes.593$total <- do.call(rbind, lapply(vcf.filtered.593, "[[", "surface"))

(top_genes.673 <- plot_top_genes(vcf.filtered.673, VI_3429_673_somatic, n = 10, "Patient 2"))
top_genes.673$total <- do.call(rbind, lapply(vcf.filtered.673, "[[", "surface"))

top.grid <- ggarrange(top_genes.593$plot, top_genes.673$plot, ncol = 3, nrow = 1, widths = c(1,1, .3),
          common.legend = T, legend = "right")

# Save the plot
ggsave(plot = top.grid, filename = file.path(date, plots_dir, "top_mutated_genes.png"),
       device = "png", width = 18, height = 6, units = "in", dpi = 300)

# 5.) Examine the expression of the shared variants
# Load expression data
readcounts <- read.csv(file = file.path("..","rnaseq",data,"readcounts.csv"),
                       sep = ",", header = T, na.strings = NA, row.names = 1)
readcounts <- readcounts[,!grepl("2DH", colnames(readcounts))]
# ------------ Patient 1 ----------------------------------------------------- #
deseq <- DESeqDataSetFromMatrix(readcounts, colData = as.data.frame(conditions),
                                design = ~1)
# Perform differential expression analysis
dds <- DESeq(deseq)
rld <- rlog(dds)

# ------------ Patient 2 ----------------------------------------------------- #
heatmap.df <- make_matrix(assay(rld), c(mapIds(org.Hs.eg.db, as.character(unique(shared_SNPs.593$geneSymbol)),
                                                keytype = "SYMBOL", column = "ENSEMBL"),
                                         mapIds(org.Hs.eg.db, as.character(unique(shared_SNPs.673$geneSymbol)),
                                                keytype = "SYMBOL", column = "ENSEMBL")))
rownames(heatmap.df) <- c(mapIds(org.Hs.eg.db, as.character(row.names(heatmap.df)),
                           keytype = "ENSEMBL", column = "SYMBOL"))

#Save image
png(file.path(date, plots_dir, "heatmap.png"), width = 10, height = 10,
    units = "in", res = 300, pointsize = 12, bg = "white")
pheatmap::pheatmap(heatmap.df, cluster_rows = T, cluster_cols = T, show_rownames = T,
                  show_colnames = T, fontsize = 8, border_color = NA, cellwidth = 10,
                  cellheight = 10)
dev.off()

# 6.) Make lollipops of common regulated proteins
# Load domain data
gff = system.file("extdata", "protein_domains.RDs", package = "maftools")
gff = readRDS(file = gff)
data.table::setDT(x = gff)

# ------------ Patient 1 ----------------------------------------------------- #
shared_SNPs.593 <- do.call(rbind, lapply(vcf.filtered.593, "[[", "surface")) %>% 
  group_by(geneSymbol, CHROM, REF, ALT, rsID, sample, POS) %>% 
  dplyr::filter(Consequence == "missense_variant") %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::group_by(geneSymbol) %>%
  dplyr::mutate(n = n()) %>% 
  dplyr::filter(n == 7) %>% 
  dplyr::select(!c(sample, n)) %>%
  dplyr::distinct(. ,.keep_all = T) %>% 
  dplyr::mutate(variant_ID = paste0(CHROM, "-", POS, "-", REF, "-", ALT))
  
  
C11orf24.data <- data.table(
  HGNC = "C11orf24",
  refseq.ID = "NM_022338",
  protein.ID = "NP_071733",
  aa.length = 449,
  Start = NA,
  End = NA,
  domain.source = NA,
  Label = NA,
  domain.anno = NA,
  pfam = NA,
  Description = NA,
  Description.1 = NA
)


C11orf24 <- data.frame(pos = c(150),
                       count = c(1),
                       Variant_Classification = "Missense_Mutation")
png(file.path(date, plots_dir, "C11orf24.png"), width = 12, height = 5,
    units = "in", res = 300, pointsize = 12, bg = "white")
make_lolliPlot(data = C11orf24, gene = "C11orf24",prot = C11orf24.data,
               refSeqID = "NM_022338", proteinID = "NP_071733")
dev.off()

CATSPER2.data <- data.table(
  HGNC = "CATSPER2",
  refseq.ID = "NM_172095",
  protein.ID = "NP_742093",
  aa.length = 530,
  Start = 172,
  End = 340,
  domain.source = "db_xref",
  Label = "Ion_trans",
  domain.anno = "Ion transport protein; pfam00520",
  pfam = "pfam00520",
  Description = "Ion transport protein",
  Description.1 = NA
)
CATSPER2 <- data.frame(pos = c(57),
                       count = c(1),
                       Variant_Classification = "Missense_Mutation")

png(file.path(date, plots_dir, "CATSPER2.png"), width = 12, height = 5,
    units = "in", res = 300, pointsize = 12, bg = "white")
make_lolliPlot(data = CATSPER2, gene = "CATSPER2", prot = CATSPER2.data,
               refSeqID = "NM_172095", proteinID = "NP_742093")
dev.off()

CATSPER2 <- data.frame(pos = c(57),
                       count = c(1),
                       Variant_Classification = "Missense_Mutation")

png(file.path(date, plots_dir, "CATSPER2.png"), width = 12, height = 5,
    units = "in", res = 300, pointsize = 12, bg = "white")
make_lolliPlot(data = CATSPER2, gene = "CATSPER2", prot = CATSPER2.data,
               refSeqID = "NM_172095", proteinID = "NP_742093")
dev.off()
# ------------ Patient 2 ----------------------------------------------------- #
shared_SNPs.673 <- do.call(rbind, lapply(vcf.filtered.673, "[[", "surface")) %>% 
  group_by(geneSymbol, CHROM, REF, ALT, rsID, sample, POS) %>% 
  dplyr::filter(Consequence == "missense_variant") %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::group_by(geneSymbol) %>%
  dplyr::mutate(n = n()) %>% 
  dplyr::filter(n == 7) %>% 
  dplyr::select(!c(sample, n)) %>%
  dplyr::distinct(. ,.keep_all = T) %>% 
  dplyr::mutate(variant_ID = paste0(CHROM, "-", POS, "-", REF, "-", ALT))

# CSMD3      chr8  G     A     "rs770488162&COSV52246013" 112859206 chr8-112859206-G-A
CSMD3 <- data.frame(pos = c(565),
                     count = c(1),
                     Variant_Classification = "Missense_Mutation")

png(file.path(date, plots_dir, "CSMD3.png"), width = 12, height = 5,
    units = "in", res = 300, pointsize = 12, bg = "white")
lollipopPlot(data = CSMD3, gene = "CSMD3",
             refSeqID = "NM_052900", proteinID = "NP_443132")
dev.off()

#DISP1      chr1  C     T     "COSV52659061"             222943251 chr1-222943251-C-T
DISP1 <- data.frame(pos = c(143),
                     count = c(1),
                     Variant_Classification = "Missense_Mutation")

png(file.path(date, plots_dir, "DISP1.png"), width = 12, height = 5,
    units = "in", res = 300, pointsize = 12, bg = "white")
lollipopPlot(data = DISP1, gene = "DISP1",
             refSeqID = "NM_032890", proteinID = "NP_116279")
dev.off()

# PCDHB3     chr5  G     A     "rs374613377&COSV50574001" 141100998 chr5-141100998-G-A
PCDHB3 <- data.frame(pos = c(117),
                    count = c(1),
                    Variant_Classification = "Missense_Mutation")

png(file.path(date, plots_dir, "PCDHB3.png"), width = 12, height = 5,
    units = "in", res = 300, pointsize = 12, bg = "white")
lollipopPlot(data = PCDHB3, gene = "PCDHB3",
             refSeqID = "NM_018937", proteinID = "NP_061760")
dev.off()

# SLC12A4    chr16 C     T     "rs781261547"               67946335 chr16-67946335-C-T
SLC12A4 <- data.frame(pos = c(815),
                     count = c(1),
                     Variant_Classification = "Missense_Mutation")

png(file.path(date, plots_dir, "SLC12A4.png"), width = 12, height = 5,
    units = "in", res = 300, pointsize = 12, bg = "white")
lollipopPlot(data = SLC12A4, gene = "SLC12A4",
             refSeqID = "NM_001145961", proteinID = "NP_001139433")
dev.off()


