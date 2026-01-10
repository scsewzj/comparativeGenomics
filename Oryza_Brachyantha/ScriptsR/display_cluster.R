#!/usr/bin/env Rscript

# ========================
# Load libraries
# ========================
library(ggplot2)
library(tools)

# ========================
# Parse command-line args
# ========================
args <- commandArgs(trailingOnly = TRUE)

# MCL cluster files
MCL_files <- if (length(args) >= 1) {
  args
} else {
  c("clusters/strict_mcl.tabular",
    "clusters/moderate_mcl.tabular")
}

plot_dir <- "plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ========================
# Function: process one MCL file
# ========================
process_mcl_file <- function(mcl_file, plot_dir) {

  clusters <- readLines(mcl_file)

  cluster_sizes <- sapply(clusters, function(line) {
    length(unlist(strsplit(line, "\\s+")))
  })

  df <- data.frame(cluster_size = cluster_sizes)

  # ---- Print summary ----
  cat("\nFile:", mcl_file, "\n")
  cat("Number of clusters:", length(cluster_sizes), "\n")
  cat("Min size:", min(cluster_sizes), "\n")
  cat("Max size:", max(cluster_sizes), "\n")
  cat("Mean size:", mean(cluster_sizes), "\n")
  cat("Median size:", median(cluster_sizes), "\n\n")

  # ---- Plot ----
  p <- ggplot(df, aes(x = cluster_size)) +
    geom_histogram(binwidth = 1, fill = "steelblue") +
    labs(
      title = paste("Cluster Size Distribution:", basename(mcl_file)),
      x = "Number of proteins per cluster",
      y = "Number of clusters"
    ) +
    theme_minimal()

  # ---- Save plot ----
  plot_name <- paste0(
    file_path_sans_ext(basename(mcl_file)),
    "_cluster_size_distribution.png"
  )

  ggsave(
    filename = file.path(plot_dir, plot_name),
    plot = p,
    width = 6,
    height = 4
  )
}

# ========================
# Run for all input files
# ========================
for (mcl_file in MCL_files) {
  process_mcl_file(mcl_file, plot_dir)
}

# =================================================================
# OPTIONAL NETWORK SECTION (kept separate on purpose)
# =================================================================
# Uncomment if needed

# library(igraph)
#
# edges_file <- "clusters/moderate_pairs.tsv"
# edges <- read.table(edges_file, header = FALSE)
# colnames(edges) <- c("protein1", "protein2", "weight")
#
# g <- graph_from_data_frame(edges, directed = FALSE)
# g <- delete_edges(g, E(g)[weight < 50])
#
# png(file.path(plot_dir, "protein_network.png"), width = 800, height = 800)
# plot(
#   g,
#   vertex.size = 5,
#   vertex.label = NA,
#   edge.width = E(g)$weight / 50,
#   edge.color = "grey",
#   main = "Protein similarity network"
# )
# dev.off()
