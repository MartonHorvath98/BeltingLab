# Set downstream path
wd <- getwd()
# Create the dated results folder
isoform_folder <- "isoform"
if (!dir.exists(file.path(data, isoform_folder))) {
  dir.create(file.path(data, isoform_folder)) # create the data folder
} 
plots_dir <- "plots" # plots directory
results_dir <- "results" # results directory
date <- format(Sys.Date(), "%Y-%m-%d") # get the current date
if (!dir.exists(file.path(date))) {
  dir.create(file.path(date)) # create the dated results folder
  dir.create(file.path(date, results_dir)) # create the results folder
  dir.create(file.path(date, plots_dir)) # create the plots folder
}

isoform.path <- file.path("..", "..", "Results", "07_isoforms")
isoform.files <- list.files(isoform.path, pattern = "^VI", full.names = T)
file.copy(from = isoform.files, recursive = T, to = file.path(data, isoform_folder))

count.files <- list.files(isoform.path, pattern = "counts", full.names = T)
file.copy(from = count.files, recursive = T, to = file.path(data, isoform_folder))

samples <- list.files(isoform.path, pattern = "^VI", full.names = F)
coldata <- data.frame("samplenames" = samples) %>%
  # extract information from the sample names
  dplyr::mutate(samplenames = as.factor(samplenames)) %>% 
  dplyr::mutate(patient = dplyr::case_when( # extract the patient IDs
    stringr::str_detect(samples, "593") ~ "593",
    stringr::str_detect(samples, "673") ~ "673"
  ),
  patient = factor(patient, levels = c("593","673"))) %>%
  dplyr::mutate(dimension = dplyr::case_when( # extract the dimension
    stringr::str_detect(samples, "-2D") ~ "2D",
    stringr::str_detect(samples, "-3D") ~ "3D",
    TRUE ~ "Tumour"
  ),
  dimension = factor(dimension, levels = c("Tumour","2D","3D"))) %>% 
  dplyr::mutate(oxygen = dplyr::case_when( # extract the growth conditions
    stringr::str_detect(samples, "H-") ~ "hypoxia",
    stringr::str_detect(samples, "N-") ~ "normoxia",
    TRUE ~ "physioxia"
  ),
  oxygen = factor(oxygen, levels = c("physioxia","normoxia","hypoxia"))) %>%
  dplyr::mutate(run = dplyr::case_when( # extract the technical replicates
    stringr::str_detect(samples, "-1-") ~ "run1",
    stringr::str_detect(samples, "-2-") ~ "run2",
    stringr::str_detect(samples, "-3-") ~ "run3"
  ),
  run = factor(run, levels = c("run1","run2","run3"))) %>%
  dplyr::mutate(
    condition = factor( # set up the factor describing the experimental design
      str_glue("{oxygen}.{dimension}"),
      levels = c("normoxia.2D","hypoxia.2D","physioxia.3D","physioxia.Tumour")))

bg = ballgown(dataDir=file.path(data, isoform_folder),
              samplePattern='VI-3429', meas='all',
              pData=coldata)

# save RData
dir.create("./RData", showWarnings = FALSE)
save(bg, file = "./RData/isoform_ballgown.RData")
rm(bg)

gene_expression = gexpr(bg)

plotTranscripts("ENSG00000146648.19", gown=bg, 
                samples=c('VI-3429-593-tumor-tissue','VI-3429-673-tumor-tissue'),
                meas='FPKM', colorby='transcript', 
                main='EGFR FPKM')

