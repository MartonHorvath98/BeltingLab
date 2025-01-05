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

# 3.) Load analysis data              
fusion.files <- list.files(file.path("../../results/06_fusion/"), pattern = "^VI.*\\.tsv", full.names = T)
fusion.files <- fusion.files[!grepl("2DH", fusion.files)]
# copy files to data directory
file.copy(fusion.files, file.path(data), overwrite = T)
# read files
fusion.files <- list.files(data, pattern = "^VI.*\\.tsv", full.names = T)
file.names <- list.files(data, pattern = "^VI.*\\.tsv", full.names = F)
file.names <- sub(pattern = ".tsv", replacement = "", file.names)
# input is a tab-delimited csv, with quantification parameters in the first row 
fusion.df <- lapply(fusion.files, read.table, stringsAsFactors=F, sep="\t", 
                    header=T, comment.char="", quote="")
names(fusion.df) <- file.names

library(dplyr)
library(tidyr)

for (i in file.names){
  fusion.df[[i]] <- fusion.df[[i]] %>% 
    dplyr::mutate(condition = i)
}


fusion.breakpoints <- lapply(fusion.df, function(x){
  return( 
    x %>%
      dplyr::select(c(1,2,3,8,10,28)) %>% 
      tidyr::separate(.,
                      col = X.FusionName, 
                      into = c("donor", "acceptor"), 
                      sep = "--") %>% 
      dplyr::mutate(FusionName = stringr::str_glue("{donor}--{acceptor}"),
                    FusionName = as.factor(FusionName)) %>% 
      tidyr::separate(., 
                      col = LeftBreakpoint, 
                      into = c("donorChr", "donorBreakpoint", "donorStrand"), 
                      sep = ":") %>% 
      tidyr::separate(., 
                      col =  RightBreakpoint, 
                      into = c("acceptorChr", "acceptorBreakpoint", "acceptorStrand"), 
                      sep = ":") %>% 
      dplyr::mutate(supportingReads = JunctionReadCount + SpanningFragCount,
                    donorBreakpoint = as.numeric(donorBreakpoint),
                    acceptorBreakpoint = as.numeric(acceptorBreakpoint)) %>% 
      dplyr::select(c(11,1,5,6,7,2,8,9,10,12,13,3,4))
  )
}) %>% do.call(rbind, .)


fusion.breakpoints <- fusion.breakpoints %>% 
  dplyr::mutate(patient = ifelse(grepl("593", condition), "593", "673"),
                patient = factor(patient, levels = c("593", "673")),
                condition = factor(condition)) %>% 
  dplyr::relocate(where(is.factor), .before = everything()) %>% 
  dplyr::group_split(patient)


names(fusion.breakpoints) <- c("593","673")

fusion.breakpoints <- lapply(fusion.breakpoints, function(x){
  x %>% 
    dplyr::group_by(FusionName) %>% 
    dplyr::mutate(start = ifelse(donorStrand == "+", 
                                 donorBreakpoint, acceptorBreakpoint),
                  start = min(start),
                  end = ifelse(donorStrand == "-", donorBreakpoint, 
                               acceptorBreakpoint),
                  end = min(end),
                  supportingReads = sum(supportingReads)) %>% 
    dplyr::select(donorChr, start, acceptorChr, end, supportingReads) %>% 
    dplyr::distinct(., .keep_all = T) %>% 
    dplyr::arrange(donorChr, start)
})



library(vcfR)
vcf.673 <- read.vcfR(file.path("A:/GBM/RStudio/Fusions/data/VI-3430-673-tumor-tissue_Manta.vcf"))

sv_data <- vcf.673@fix %>% 
  as.data.frame() %>%
  dplyr::select(CHROM, POS, INFO)
sv_data <- sv_data %>%
  mutate(
    SVTYPE = gsub(".*SVTYPE=([A-Z]+);.*", "\\1", INFO),
    END = as.numeric(gsub(".*END=([0-9]+);.*", "\\1", INFO)),
    CHROM = as.character(CHROM)
  )

circos_data <- sv_data %>%
  filter(!is.na(SVTYPE) & !is.na(END) & 
           CHROM %in% paste0("chr", c(1,2,3,5,6,7,9,12,13,14,19,"X"))) %>%
  dplyr::select(CHROM, POS, END, SVTYPE) %>% 
  dplyr::mutate(CHROM = as.factor(CHROM),
                POS = as.numeric(POS),
                END = as.numeric(END),
                SVTYPE = as.factor(SVTYPE))


library(circlize)
svg(file.path(date, plots_dir, "circos_SV_593.svg"), width = 10, height = 10)
circos.initializeWithIdeogram(plotType = c("axis", "labels"))
text(0, 0, "Patient 1", cex = 2)
circos.genomicLink(
  region1 = circos_data[, c("CHROM", "POS", "POS")],
  region2 = circos_data[, c("CHROM", "END", "END")],
  col = ifelse(circos_data$SVTYPE == "DUP", "steelblue", ifelse(circos_data$SVTYPE == "DEL", "violet", "salmon")),
  border = NA, h = 1
)
circos.genomicLink(
  region1 = fusion.breakpoints[["593"]][,c("donorChr", "start", "start")],
  region2 = fusion.breakpoints[["593"]][,c("acceptorChr", "end", "end")],
  col = "red", lw = 4, h = .7
)
circos.labels(sectors = fusion.breakpoints[["593"]]$donorChr,
              x = fusion.breakpoints[["593"]]$start, niceFacing = T,
              labels = strwrap(fusion.breakpoints[["593"]]$FusionName, width = 20),
              side = "inside", cex = 1, col = "black")
dev.off()
