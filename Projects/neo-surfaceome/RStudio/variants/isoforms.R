# Set downstream path
wd <- getwd()
folder <- "mRNA"
# Create the dated results folder
data_dir <- "data"
if (!dir.exists(file.path(data_dir, folder))) {
  dir.create(file.path(data_dir, folder)) # create the main results folder
}
results_dir <- "results"
if (!dir.exists(file.path(results_dir, folder))) {
  dir.create(file.path(results_dir, folder)) # create the main results folder
}
# get the current date
date <- format(Sys.Date(), "%Y-%m-%d")
#plots directory
plots_dir <- "plots"
#tables directory
tables_dir <- "tables"
if (!dir.exists(file.path(results_dir, folder, date))) {
  dir.create(file.path(results_dir, folder, date), recursive = T) # create the dated results folder
  dir.create(file.path(results_dir, folder, date, tables_dir)) # create the tables folder
  dir.create(file.path(results_dir, folder, date, plots_dir)) # create the plots folder
}


mrna.files <- list.files(mrna.path, pattern = "counts.txt$", full.names = T)
file.copy(from = mrna.files, to = file.path(data_dir, folder))



isoformQuant <- importIsoformExpression(
  parentDir = file.path("..", "..", "Results", "07_isoform")
)