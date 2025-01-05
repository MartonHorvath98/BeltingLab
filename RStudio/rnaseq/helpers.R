library(msigdbr)
# load databases
msigdbr_df <- msigdbr(species = "Homo sapiens") # save local databases
# extract KEGG- and REACTOME pathways from MSigDB
pathways <- msigdbr_df %>% 
  dplyr::filter(
    gs_cat == "C2", # only canonical representations (compiled by experts)
    gs_subcat %in% c("CP:KEGG", "CP:REACTOME") # KEGG and Reactome pathways
  )
# extract GO terms and Hallmark genesets
terms <- msigdbr_df %>% 
  dplyr::filter(
    gs_cat == "C5" & # Ontology gene sets
      gs_subcat == "GO:BP" | # GO terms, excluding phenotype ontology
      gs_cat == "H" # Hallmark gene sets
)

# Prepare the SURFME filter
library(openxlsx)
surfme <- read.xlsx(file.path(data,"surfme_v2023.xlsx"))
surfme_categories <- c("GPI","SinglePass","MultiPass","Cellular_Membrane",
                       "Extracellular_Domain","GOplasmaMb","GOCellSurf",
                       "GOexternal","Surfy")
surfme <- surfme %>% 
  tidyr::separate_rows(Gene.names1, sep = " ") %>% 
  dplyr::rename("UniprotID" = Entry, 
                "Description" = Protein.names1,
                "geneID" = Gene.names1) %>% 
  dplyr::mutate(EnsemblID = mapIds(org.Hs.eg.db, keys = UniprotID, column = "ENSEMBL", 
                                   keytype = "UNIPROT", multiVals = "first")) %>%
  dplyr::select(c(1,2,4,15, 6:10,14,11,12,13)) %>% 
  dplyr::mutate(across(all_of(surfme_categories), ~ ifelse(is.na(.), F, T))) %>% 
  dplyr::distinct(EnsemblID, .keep_all = T) %>%
  dplyr::filter(complete.cases(.)) 

rm("surfme_categories")