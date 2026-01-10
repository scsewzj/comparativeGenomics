#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

# ========================
# Parse command-line args
# ========================
args <- commandArgs(trailingOnly = TRUE)

# if (length(args) < 1) {
#   stop("Usage: KaKs_all_ks_input_<stringency>.tsv <input.tsv> [output.pdf]")
# }

input_file  <- ifelse(length(args) >= 1, args[1], "./KaKs_all_ks_input_moderate.tsv")
output_file <- ifelse(length(args) >= 2, args[2], "dS_evolution_moderate.png")

plot_dir <- "plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)


# Read the data (example: tab-delimited file)
# pair    dS
# OB01G10050.1_OB01G34280.1  1.1700
# OB01G10050.1_OB05G10330.1  1.2456
# OB01G10070.1_OB01G16230.1  2.0768

data <- read.table(input_file, header = TRUE, sep = "\t")
data <- data[is.finite(data$dS) & is.finite(data$dS_SE) & data$dS > 0 & data$dS < 6 & data$dS_SE < 5,c("pair","dS")]


# Define bin width
binwidth <- 0.01

# Create breaks
breaks <- seq(
  from = floor(min(data$dS) / binwidth) * binwidth,
  to   = ceiling(max(data$dS) / binwidth) * binwidth,
  by   = binwidth
)

# Plot histogram

p <- ggplot(data, aes(x = data$dS)) +
  geom_histogram(
    breaks = breaks,
    fill = "steelblue", # steelblue # darkgreen 
  ) +
  labs(
    title = "(L) O. brachyantha",
    x = "dS",
    y = "Number of gene pairs"
  ) +
  theme_minimal()

plot(p)

ggsave(file.path(plot_dir, output_file), plot = p, width = 7, height = 5)

message("Plot saved to: ", file.path(plot_dir, output_file))

# # install.packages("pracma")  # if not installed
# library(pracma)
# 
# # Apply to histogram counts
# hist_counts <- hist(data$dS, breaks = breaks, plot = FALSE)
# peaks <- findpeaks(hist_counts$counts, nups = 1, ndowns = 1)  # simple peak detection
# print(peaks)  # gives height, location index, etc.
