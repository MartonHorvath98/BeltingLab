################################################################################
# Created: 2024 05 21 ; Last Modified: 2025 06 23 ; KGO & MH                   #
################################################################################
#                              Data Processing                                 #
################################################################################
# Set working directory
wd <- getwd()
# Load packages
if (file.exists(file.path(wd, "packages.R"))) {
  source(file = file.path(wd, "packages.R"))
} else {
  stop("Required file 'packages.R' not found in the working directory.")
}
# Source data processing functions
if (file.exists(file.path(wd, "functions.R"))) {
  source(file = file.path(wd, "functions.R"))
} else {
  stop("Required file 'functions.R' not found in the working directory.")
}

# Export session info and save to file
sink(file = "../sessionInfo.txt")
devtools::session_info()
sink()

# ---------------------------------------------------------------------------- #
# -   Svenja's CC +LD and CC no LD (Clariom D Human Pico) Affymetrix         - #
# ---------------------------------------------------------------------------- #
data_dir <- "../data/raw/Affymetrix/CCLD/"
CCLD.dat <- read.celfiles(list.celfiles(data_dir, full.names = TRUE))

# Transcript filtering
CCLD.expr <- normalizeTranscript(CCLD.dat, clariomdhumantranscriptcluster.db)
# log2Foldchange--> log2(exp)-log2(ctr)
CCLD.expr$log2FoldChange <- CCLD.expr$`CC_LD_(Clariom_D_Human).CEL` - CCLD.expr$`CC_noLD_(Clariom_D_Human).CEL`

####### Remove duplicate symbols based on lgFC
CCLD.expr <- removeDuplicates(.data = CCLD.expr,.column = "log2FoldChange",
                              .symbol = "SYMBOL")

####### Check how many up and down- regulated genes
length(which(CCLD.expr$log2FoldChange >=0.5)) # 6082
length(which(CCLD.expr$log2FoldChange <=(-0.5))) # 897

# save RData
dir.create("./RData", showWarnings = FALSE)
save(CCLD.expr, file = "./RData/CCLD_processedData.RData")
# save file
dir.create("../data/processed/expression_matrices/", showWarnings = FALSE)
write.xlsx(CCLD.expr, file = "../data/processed/expression_matrices/CCLD_processedData.xlsx")

# ---------------------------------------------------------------------------- #
# -      Hugo's Primary cells 2D vs 3D (Clariom D Human Pico) Affymetrix     - #
# ---------------------------------------------------------------------------- #
data_dir <- "../data/raw/Affymetrix/HGCC/"
HGCC.dat <- read.celfiles(list.celfiles(data_dir, full.names = TRUE))

# Transcript filtering
HGCC.expr <- normalizeTranscript(HGCC.dat, clariomdhumantranscriptcluster.db)

# Compare cell lines
HGCC.deg <- limmaDEA(.data = HGCC.expr,
                     .design = c("U3017_2D", "U3017_2D", "U3017_2D",
                                 "U3017_3D", "U3017_3D", "U3017_3D",
                                 "U3047_2D", "U3047_2D", "U3047_2D",
                                 "U3047_3D", "U3047_3D", "U3047_3D",
                                 "U3054_2D", "U3054_2D", "U3054_2D",
                                 "U3054_3D", "U3054_3D", "U3054_3D"),
                     .contrast = c("U3017_3D-U3017_2D",
                                   "U3047_3D-U3047_2D",
                                   "U3054_3D-U3054_2D")) 

names(HGCC.deg) <- c("U3017", "U3047", "U3054")
# Compare global
HGCC.deg <- c(HGCC.deg,
              "global" = limmaDEA(.data = HGCC.expr,
                                  .design = c("2D","2D","2D","3D","3D","3D",
                                              "2D","2D","2D","3D","3D","3D",
                                              "2D","2D","2D","3D","3D","3D"),
                                  .contrast = c("X3D - X2D")))

# Remove duplicate IDs
HGCC.deg <- lapply(HGCC.deg, removeDuplicates, .column = "t", .symbol = "ID.Symbol")

HGCC.df <- merge.rec(HGCC.deg, by = c("ID.ID"), all = T, suffix = c("",""))
HGCC.df <- HGCC.df[,c(1:3, 5, # Identifier
                      4, 6:8, # U3017 2D vs 3D 
                      12, 14:16, # U3047 2D vs 3D 
                      20, 22:24, # U3054 2D vs 3D 
                      28, 30:32)] %>%  # Global 2D vs 3D
           setNames(., c("PROBEID", "Symbol", "EntrezID", "AveExpr",
             "log2FoldChange.U3017", "statT.U3017", "pvalue.U3017", "padj.U3017",
             "log2FoldChange.U3047", "statT.U3047", "pvalue.U3047", "padj.U3047",
             "log2FoldChange.U3054", "statT.U3054", "pvalue.U3054", "padj.U3054",
             "log2FoldChange.global", "statT.U3054", "pvalue.global", "padj.global"))
                  

####### Check how many up and down- regulated genes
length(which(HGCC.df$log2FoldChange.global >= 0.5)) # 2539
length(which(HGCC.df$log2FoldChange.global <= (-0.5))) # 1831

# Save RData
save(HGCC.expr, HGCC.deg, file = "./RData/HGCC_processedData.RData")
# Save file
sapply(names(HGCC.deg), function(x){
  write.xlsx(HGCC.deg[[x]], file = paste0("../data/processed/expression_matrices/HGCC_",x,"_processedData.xlsx"))
})

# ---------------------------------------------------------------------------- #
# -     U87 Chronic Acidosis AA vs NA & HOX vs NOX (Illumina BeadChip)       - #
# ---------------------------------------------------------------------------- #
data_dir <- "../data/raw/Illumina/U87/"
U87.dat <- read.ilmn(files=file.path(data_dir, "Sample_Probe_Summary.txt"),
                     ctrlfiles=file.path(data_dir, "Control_Probe_Summary.txt"))

samples <- c(
  # Chronic acidosis
  "200118400068_I" = "control_sel", #13U87selctrlpH74 - Selection control pH 7.4
  "200118400035_I" = "control_sel", #14U87selctrlpH74 - Selection control pH 7.4
  "200118400035_D" = "control_sel", #15U87selctrlpH74 - Selection control pH 7.4
  "200118400035_K" = "sel_pH647", #10U87selpH647 - Selection pH 6.4
  "200118400035_G" = "sel_pH647", #11U87selpH647 - Selection pH 6.4
  "200118400035_F" = "sel_pH647", #12U87selpH647 - Selection pH 6.4
  # Acute acidosis
  "200118400068_K" = "control_acu", #7U87aactrl74 - Acute acidosis control pH 7.4
  "200118400068_D" = "control_acu", #8U87aactrl74 - Acute acidosis control pH 7.4
  "200118400068_A" = "control_acu", #9U87aactrl74 - Acute acidosis control pH 7.4
  "200118400068_C" = "acu_pH68", #4U87aapH68 - Acute acidosis pH 6.8
  "200118400068_G" = "acu_pH68", #5U87aapH68 - Acute acidosis pH 6.8
  "200118400033_E" = "acu_pH68", #6U87aapH68 - Acute acidosis pH 6.8
  "200118400033_B" = "acu_pH64", #1U87aapH64 - Acute acidosis pH 6.4
  "200118400068_B" = "acu_pH64", #2U87aapH64 - Acute acidosis pH 6.4
  "200118400035_B" = "acu_pH64", #3U87aapH64 - Acute acidosis pH 6.4
  # Normoxia
  "200118400033_G" = "control_nox", #nox3 - atmospheric O2 (21%) control
  "200118400035_L" = "control_nox", #nox2 - atmospheric O2 (21%) control
  "200118400033_H" = "control_nox", #nox1 - atmospheric O2 (21%) control
  # Hypoxia
  "200118400035_E" = "hypoxia", #hox3 - 1% O2
  "200118400035_H" = "hypoxia", #hox2 - 1% O2
  "200118400033_I" = "hypoxia") #hox1 - 1% O2

# Connect to DB
dbCon=org.Hs.eg_dbconn()

# SQL query
sqlQuery='SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'

# Query the database
aliasSymbol=dbGetQuery(dbCon, sqlQuery)

#normalization and background correction
U87.expr <- normalizeIllumina(U87.dat, .dbconn = aliasSymbol, .samples = names(samples))

#Perform DEG analysis
U87.deg <- limmaDEA(.data = U87.expr,
                    .design = samples,
                    .contrast = c("sel_pH647-control_sel",
                                  "acu_pH68-control_acu",
                                  "acu_pH64-control_acu",
                                  "hypoxia-control_nox"))

# Remove duplicate IDs
names(U87.deg) <- c("sel_pH647-control_sel", "acu_pH68-control_acu",
                    "acu_pH64-control_acu", "hypoxia-control_nox")
U87.deg <- lapply(U87.deg, removeDuplicates, .column = "t", .symbol = "ID.Symbol")

####### Check how many up and down- regulated genes
length(which(U87.deg$`sel_pH647-control_sel`$logFC >= 0.5)) # 1876
length(which(U87.deg$`sel_pH647-control_sel`$logFC <= (-0.5))) # 1866

# Save RData
save(U87.expr, U87.deg, file = "./RData/U87_processedData.RData")
# Save file
sapply(names(U87.deg), function(x){
  write.xlsx(U87.deg[[x]], file = paste0("../data/processed/expression_matrices/U87_",x,"_processedData.xlsx"))
})

# ---------------------------------------------------------------------------- #
# -   PANC1 Chronic Acidosis AA vs NA (Clariom D Human Pico) Affymetrix      - #
# ---------------------------------------------------------------------------- #
data_dir <- "../data/raw/Affymetrix/PANC1/"
PANC1.dat <- read.celfiles(list.celfiles(data_dir, full.names = TRUE))

# Transcript filtering
PANC1.expr <- normalizeTranscript(PANC1.dat, clariomdhumantranscriptcluster.db)

# Perform DEG analysis
PANC1.deg <- limmaDEA(.data = PANC1.expr,
                     .design = c("PANC1_NA", "PANC1_NA","PANC1_NA",
                                 "PANC1_AA", "PANC1_AA", "PANC1_AA"),
                     .contrast = c("PANC1_AA - PANC1_NA"))

PANC1.deg <- removeDuplicates(PANC1.deg[[1]], .column = "t", .symbol = "ID.Symbol")

####### Check how many up and down- regulated genes
length(which(PANC1.deg$logFC >= 0.5)) # 1038
length(which(PANC1.deg$logFC <= (-0.5))) # 659

# Save RData
save(PANC1.expr, PANC1.deg, file = "./RData/PANC1_processedData.RData")
# Save file
write.xlsx(PANC1.deg, file = paste0("../data/processed/expression_matrices/PANC1_processedData.xlsx"))


################################################################################
# Clean up the environment                                                     #
################################################################################
# Close any open connections
if (exists("dbCon") && dbIsValid(dbCon)) {
  dbDisconnect(dbCon)
}
# Run garbage collection to free up memory
gc() 
# Clear the environment
rm(list = ls()) 




