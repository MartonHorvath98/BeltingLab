conditions <- c("VI-3429-593-2DN-1","VI-3429-593-2DN-2","VI-3429-593-2DN-3","VI-3429-593-3D-1",
                "VI-3429-593-3D-2","VI-3429-593-3D-3","VI-3429-593-tumor-tissue", "VI-3429-673-2DN-1",
                "VI-3429-673-2DN-2","VI-3429-673-2DN-3","VI-3429-673-3D-1","VI-3429-673-3D-2",
                "VI-3429-673-3D-3","VI-3429-673-tumor-tissue")


coldata <- data.frame("samplenames" = conditions) %>%
  # extract information from the sample names
  dplyr::mutate(samplenames = as.factor(samplenames)) %>% 
  dplyr::mutate(patient = dplyr::case_when( # extract the patient IDs
    stringr::str_detect(conditions, "593") ~ "593",
    stringr::str_detect(conditions, "673") ~ "673"
  ),
  patient = factor(patient, levels = c("593","673"))) %>%
  dplyr::mutate(dimension = dplyr::case_when( # extract the dimension
    stringr::str_detect(conditions, ".2D") ~ "2D",
    stringr::str_detect(conditions, ".3D") ~ "3D",
    TRUE ~ "Tumour"
  ),
  dimension = factor(dimension, levels = c("Tumour","2D","3D"))) %>% 
  dplyr::mutate(oxygen = dplyr::case_when( # extract the growth conditions
    stringr::str_detect(conditions, "H.") ~ "hypoxia",
    stringr::str_detect(conditions, "N.") ~ "normoxia",
    TRUE ~ "physioxia"
  ),
  oxygen = factor(oxygen, levels = c("physioxia","normoxia","hypoxia"))) %>%
  dplyr::mutate(run = dplyr::case_when( # extract the technical replicates
    stringr::str_detect(conditions, ".1$") ~ "run1",
    stringr::str_detect(conditions, ".2$") ~ "run2",
    stringr::str_detect(conditions, ".3$") ~ "run3"
  ),
  run = factor(run, levels = c("run1","run2","run3"))) %>%
  dplyr::mutate(
    condition = factor( # set up the factor describing the experimental design
      str_glue("{oxygen}.{dimension}"),
      levels = c("normoxia.2D","hypoxia.2D","physioxia.3D","physioxia.Tumour")))

coldata <- coldata %>% 
  dplyr::mutate(condition = relevel(condition, ref = "normoxia.2D")) %>% 
  split.data.frame(.$patient, drop = T)

# Define the conversion matrix
conv = c("T>C", "T>C", "C>T", "C>T", "T>A", "T>A", "T>G", 
         "T>G", "C>A", "C>A", "C>G", "C>G")
names(conv) = c("A>G", "T>C", "C>T", "G>A", "A>T", "T>A", 
                "A>C", "T>G", "C>A", "G>T", "C>G", "G>C")
conv.class = c("Ti", "Ti", "Tv", "Tv", "Tv", "Tv")
names(conv.class) = c("T>C", "C>T", "T>A", "T>G", "C>A", 
                      "C>G")

# Define a color scale for mutations
mut_colors = c("T>C" = "red", "C>T" = "orange", "T>A" = "violet", 
               "T>G" = "purple", "C>A" = "turquoise", "C>G" = "blue")

consequence_colors = c("3_prime_UTR_variant" = "#F44336","5_prime_UTR_variant" = "orange","intron_variant" = "#9C27B0","missense_variant" = "#673AB7", 
  "splice_donor_region_variant&intron_variant"="#3F51B5","synonymous_variant"= "#03A9F4")
