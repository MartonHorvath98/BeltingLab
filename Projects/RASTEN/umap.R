# Load packages
if (!requireNamespace("openxlsx", quietly = TRUE))
  install.packages("openxlsx")
suppressPackageStartupMessages(library(openxlsx))
if (!requireNamespace("umap", quietly = TRUE))
  install.packages("umap")
suppressPackageStartupMessages(library(umap))
if (!requireNamespace("igraph", quietly = TRUE))
  install.packages("igraph")
suppressPackageStartupMessages(library(igraph))

# Load and inspect data
data_path <- "data"
rasten_inputed <- read.xlsx(file.path(data_path, "RASTEN_impLCMD_Pat_log2_Normalized.xlsx"), 
                            sheet = 1, rowNames = T)

# For the sake of the comparability between baseline and cycle 3, I remove patients
# without a cycle 3 measurement point (for now)
filter <- gsub(pattern = "_B$|_C3$", replacement = "", colnames(rasten_inputed)) |> table()
filter <- names(filter[filter == 2])

rasten_inputed <- rasten_inputed[, grepl(paste(filter, collapse = "|"), colnames(rasten_inputed))]
rasten_inputed <- t(rasten_inputed)

# Split the data into baseline and cycle 3
rasten_baseline <- rasten_inputed[grepl("_B$", rownames(rasten_inputed)), ]
rasten_cycle3 <- rasten_inputed[grepl("_C3$", rownames(rasten_inputed)), ]

baseline_umap <- umap(rasten_baseline)
cycle3_umap <- umap(rasten_cycle3)

# Calculate distances between points
baseline_distances <- as.matrix(dist(baseline_umap$layout))
cycle3_distances <- as.matrix(dist(cycle3_umap$layout))

# Calculate adjacency matrices
baseline_adjacency <- baseline_distances < 0.5
cycle3_adjacency <- cycle3_distances < 0.5

# Create igraph objects
baseline_graph <- graph_from_adjacency_matrix(baseline_adjacency, mode = "undirected", diag = TRUE)
cycle3_graph <- graph_from_adjacency_matrix(cycle3_adjacency, mode = "undirected", diag = TRUE)

# Cluster the graphs
baseline_clusters <- cluster_fast_greedy(baseline_graph) |> membership()
cycle3_clusters <- cluster_fast_greedy(cycle3_graph) |> membership()

# PLot the graphs
(baseline_plot <- ggplot(setNames(as.data.frame(baseline_umap$layout), c("x","y"))) +
  geom_point(aes(x = x, y = y, color = as.factor(baseline_clusters)),
             size = 3) + 
  stat_ellipse(geom = "polygon", 
               aes(x = x, y = y, color =  as.factor(baseline_clusters)), fill = NA, 
               type = "norm", level = 0.68, linewidth = 1,
               alpha = .3, linetype = 2, show.legend = F) + 
  scale_color_manual(name = "Cluster",
                     values = viridis::magma(length(unique(baseline_clusters)), begin = 0, end = 0.85)))

(cycle3_plot <- ggplot(setNames(as.data.frame(cycle3_umap$layout), c("x","y"))) +
  geom_point(aes(x = x, y = y, color = as.factor(cycle3_clusters)),
             size = 3) + 
  stat_ellipse(geom = "polygon", 
               aes(x = x, y = y, color =  as.factor(cycle3_clusters)), fill = NA, 
               type = "norm", level = 0.68, linewidth = 1,
               alpha = .3, linetype = 2, show.legend = F) + 
  scale_color_manual(name = "Cluster",
                     values = viridis::magma(length(unique(cycle3_clusters)), begin = 0, end = 0.85)))
