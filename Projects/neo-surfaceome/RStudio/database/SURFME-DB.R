# Install and load necessary packages
source("packages.R")

### 1.)
# Create a connection to a local SQLite database file
SFcon <- dbConnect(RSQLite::SQLite(), "SURFME.sqlite")

### 2.)
# Prepare gene identifier list
SURFME <- read.table(file.path("..","..","config","surfme_genes.txt"), header = TRUE, sep = "\t")
# Add annotations
symbols <- mapIds(org.Hs.eg.db, keys = SURFME$Ensembl.gene.ID, keytype = "ENSEMBL", column = "SYMBOL")
entrez <- mapIds(org.Hs.eg.db, keys = SURFME$Ensembl.gene.ID, keytype = "ENSEMBL", column = "ENTREZID")
# Create the SURFME genes table
dbExecute(SFcon, "CREATE TABLE genes (
  gene_id INTEGER PRIMARY KEY AUTOINCREMENT,
  gene_name TEXT UNIQUE,
  gene_symbol TEXT,
  entrez_id TEXT
)")

for (gene in SURFME$Ensembl.gene.ID) {
  dbExecute(SFcon, "INSERT OR IGNORE INTO genes (gene_name, gene_symbol, entrez_id) VALUES (?, ?, ?)", 
            params = list(gene, symbols[gene], entrez[gene]))
}

### 3.)
# Add STRING identifiers
string_db <- rba_string_map_ids(SURFME$Ensembl.gene.ID, species = 9606)
# Create the SURFME genes table
dbExecute(SFcon, "CREATE TABLE string_db (
  _id INTEGER PRIMARY KEY AUTOINCREMENT,
  gene_id INTEGER,
  gene_symbol TEXT UNIQUE,
  string_id TEXT,
  description TEXT,
  FOREIGN KEY(gene_id) REFERENCES genes(gene_id)
)")

# Insert identifiers to the table
for (i in 1:nrow(string_db)) {
  gene_id <- dbGetQuery(SFcon, "SELECT gene_id FROM genes WHERE gene_name = ?", params = list(string_db$queryItem[i]))$gene_id
  dbExecute(SFcon, "INSERT INTO string_db (gene_id, gene_symbol, string_id, description) VALUES (?, ?, ?, ?)", 
            params = list(gene_id, string_db$preferredName[i], string_db$stringId[i], string_db$annotation[i]))
}

### 4.)
# STRING interaction network
string_id <- string_db$stringId
int_partners <- lapply(string_id, function(x) {
  rba_string_interaction_partners(x, species = 9606, 
                                  required_score = 900, 
                                  network_type = "physical") 
})
# Create the STRING interactions table
int_network <- do.call(rbind, int_partners)
rm(int_partners)

dbExecute(SFcon, "CREATE TABLE string_interactions (
  _id INTEGER PRIMARY KEY AUTOINCREMENT,
  gene_symbolA INTEGER,
  gene_symbolB INTEGER,
  string_idA TEXT,
  string_idB TEXT,
  score REAL,
  score_details TEXT
)")

for (i in 1:nrow(int_network)) {
  score_details <- paste0("(", "nscore: ", int_network$nscore[i], ", fscore: ", int_network$fscore[i],
                          ", pscore: ", int_network$pscore[i], ", ascore: ", int_network$ascore[i], 
                          ", escore: ", int_network$escore[i], ", dscore: ", int_network$dscore[i], 
                          ", tscore: ", int_network$tscore[i], ")") 
  dbExecute(SFcon, "INSERT INTO string_interactions (gene_symbolA, gene_symbolB, string_idA, string_idB, score, score_details) VALUES (?, ?, ?, ?, ?, ?)", 
            params = list(int_network$preferredName_A[i], int_network$preferredName_B[i],
                          int_network$stringId_A[i], int_network$stringId_B[i],
                          int_network$score[i], score_details))
}

### 5.)
# Add the STRING interactions to the database
library(biomaRt)
ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
pfam_data <- lapply(SURFME$ENSEMBL, function(x) {
  biomaRt::getBM(
    attributes = c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol", 
                   "pfam", "pfam_start", "pfam_end", # PFAM domains
                   "signalp", "signalp_start","signalp_end", # cleavage sites
                   "tmhmm", "tmhmm_start", "tmhmm_end" # transmembrane domains
                   ),
  filters = "ensembl_gene_id",
  values = "ENSG00000074621",
  mart = ensembl
)})

PFAMdb = PFAM_dbconn()
dbGetQuery(PFAMdb, "SELECT * FROM pdb")
