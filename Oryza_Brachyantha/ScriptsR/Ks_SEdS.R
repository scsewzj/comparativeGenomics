#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

# ========================
# Parse command-line args
# ========================
args <- commandArgs(trailingOnly = TRUE)

# if (length(args) < 1) {
#   stop("Usage: plot_dS_SE_vs_dS.R <input.tsv> [output.pdf]")
# }

input_file  <- ifelse(length(args) >= 1, args[1], "./KaKs_all_ks_input_moderate.tsv")
output_file <- ifelse(length(args) >= 2, args[2], "dS_SE_vs_dS.png")

plot_dir <- "plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ========================
# Load data
# ========================
df <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Basic sanity check
required_cols <- c("dS", "dS_SE")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# Remove invalid values
threshold_SE <- 100
threshold <- 6
df <- df[is.finite(df$dS) & is.finite(df$dS_SE) & df$dS >= 0 & df$dS < threshold & df$dS_SE <= threshold_SE, ]

# ========================
# Plot
# ========================
p <- ggplot(df, aes(x = dS, y = dS_SE)) +
  geom_point(alpha = 0.5, size = 1) +
  labs(
    title = paste("dS_SE as a Function of dS<",threshold),
    x = "dS (synonymous substitutions per site)",
    y = "SE(dS)"
  ) +
  theme_minimal(base_size = 12)
plot(p)
# ========================
# Save
# ========================
ggsave(file.path(plot_dir, paste0("Moderate_cutoff_", threshold, "_SE_", threshold_SE, "_", output_file)), plot = p, width = 7, height = 5)

message("Plot saved to: ", file.path(plot_dir, paste0("cutoff_", threshold, "_SE_", threshold_SE, "_", output_file)))
