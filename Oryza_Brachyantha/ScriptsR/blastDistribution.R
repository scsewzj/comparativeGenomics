#!/usr/bin/env Rscript

# Load libraries
library(ggplot2)
library(reshape2)

# ========================
# Parse command-line args
# ========================
args <- commandArgs(trailingOnly = TRUE)

blast_file  <- ifelse(length(args) >= 1, args[1], "append_all_vs_all.tsv") # Input file
header      <- ifelse(length(args) >= 2, as.logical(args[2]), FALSE)
sep_char    <- ifelse(length(args) >= 3, args[3], "")

plot_dir <- "plots"

# ------------------------
# Create output directory
# ------------------------
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# Read the BLAST-like file
blast <- read.table(blast_file, header = header, sep = sep_char, stringsAsFactors = FALSE)

# Assign column names (14 fields)
colnames(blast) <- c("Query", "Subject", "PercIdentity", "AlignLength", "Mismatches",
                     "GapOpens", "Qstart", "Qend", "Sstart", "Send", "Evalue",
                     "BitScore", "QueryLength", "SubjectLength")

# ------------------------
# Calculate coverage
# ------------------------
blast$QueryCoverage <- (blast$Qend - blast$Qstart + 1) / blast$QueryLength * 100
blast$SubjectCoverage <- (blast$Send - blast$Sstart + 1) / blast$SubjectLength * 100

# ------------------------
# Select columns to plot
# ------------------------
df_plot <- blast[, c("PercIdentity", "QueryCoverage", "SubjectCoverage")]

# Melt for ggplot
df_melt <- melt(df_plot, variable.name = "Metric", value.name = "Value")

# ------------------------
# Plot distributions with vertical line legend
# ------------------------
vline_data <- data.frame(xintercept = c(50, 70),
                         LineType = c("50%", "70%"))

p <- ggplot(df_melt, aes(x = Value, fill = Metric)) +
  geom_histogram(alpha = 0.6, bins = 30, position = "identity", color = NA) +
  facet_wrap(~Metric, scales = "free") +
  theme_minimal() +
  labs(title = "BLAST Hit Metrics Distribution",
       x = "Percentage",
       y = "Count") +
  geom_vline(data = vline_data, aes(xintercept = xintercept, color = LineType),
             linetype = "dashed", linewidth = 0.8, show.legend = TRUE) +
  scale_color_manual(name = "Thresholds", values = c("50%" = "red", "70%" = "blue")) +
  theme(legend.position = "right")

# ------------------------
# Save plot
# ------------------------
out_file <- file.path(plot_dir, "blast_distribution.png")
ggsave(out_file, p, width = 10, height = 5)

# Print plot
print(p)
