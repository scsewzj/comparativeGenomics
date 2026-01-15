#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# ========================
# Parse command-line args
# ========================
args <- commandArgs(trailingOnly = TRUE)

blast_raw   <- ifelse(length(args) >= 1, args[1], "../append_all_vs_all.tsv")
blast_filt1 <- ifelse(length(args) >= 2, args[2], "../clusters/moderate.tsv")
blast_filt2 <- ifelse(length(args) >= 3, args[3], "../clusters/strict.tsv")
header      <- ifelse(length(args) >= 4, as.logical(args[4]), FALSE)
sep_char    <- ifelse(length(args) >= 5, args[5], "")

plot_dir <- "../plots"

# Create output directory if it doesn't exist
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ========================
# Function to get unique accession IDs
# ========================
get_unique_accessions <- function(blast_file, header = FALSE, sep = "") {
  df <- read.table(blast_file, header = header, sep = sep, stringsAsFactors = FALSE)
  unique(c(df$V1, df$V2))
}

# Get unique accession IDs
raw_accessions  <- get_unique_accessions(blast_raw, header, sep_char)
filt1_accessions <- get_unique_accessions(blast_filt1, header, sep_char)
filt2_accessions <- get_unique_accessions(blast_filt2, header, sep_char)

# Define files and labels
blast_files <- list(filt1_accessions = filt1_accessions,
                    filt2_accessions = filt2_accessions)
file_labels <- c("L", "H")

# Initialize counts dataframe
counts_df <- data.frame(
  File = character(),
  Hits = integer(),
  Sequences = integer(),
  Singletons = integer(),
  stringsAsFactors = FALSE
)

# Calculate counts for each filtered file
for (i in seq_along(blast_files)) {
  # Corresponding filtered file for Hits
  blast_table <- if (i == 1) read.table(blast_filt1, header = header, sep = sep_char, stringsAsFactors = FALSE) else
    read.table(blast_filt2, header = header, sep = sep_char, stringsAsFactors = FALSE)
  
  hits <- nrow(blast_table)
  sequences <- length(blast_files[[i]])
  singletons <- sum(!blast_files[[i]] %in% raw_accessions) # Not strictly needed, should be zero
  # Actually, singletons = accessions in raw but not in filtered
  singletons <- sum(!raw_accessions %in% blast_files[[i]])
  
  counts_df <- rbind(counts_df,
                     data.frame(File = file_labels[i],
                                Hits = hits,
                                Sequences = sequences,
                                Singletons = singletons))
}

# Transform data for ggplot
counts_long <- counts_df %>%
  pivot_longer(cols = c(Hits, Sequences, Singletons),
               names_to = "Type",
               values_to = "Count")

# Plot
p <- ggplot(counts_long, aes(x = File, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", color = NA) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5) +
  scale_fill_manual(values = c("Hits" = "#6caed6", "Sequences" = "#cbaa7a", "Singletons" = "#bebebe")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(y = "# Hits", x = "", fill = "Type") +
  theme_minimal(base_size = 14) +
  ggtitle("O. brachyantha") +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(face="italic"))
plot(p)
# Save plot
ggsave(filename = file.path(plot_dir, "blast_counts.png"), plot = p, width = 6, height = 4)
