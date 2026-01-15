#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ========================
# Parse command-line args
# ========================
args <- commandArgs(trailingOnly = TRUE)
input_file   <- ifelse(length(args) >= 1, args[1], "../KaKs_all_ks_input_strict.tsv")
tags_file    <- ifelse(length(args) >= 2, args[2], "../TAGS/TAGs_0_spacers_IDs.txt")
non_tags_file<- ifelse(length(args) >= 3, args[3], "../TAGS/non_TAGs.txt")
singletons_file <- ifelse(length(args) >= 4, args[4], "../TAGS/strict.singleton_IDs.txt")
output_file  <- ifelse(length(args) >= 5, args[5], "dS_evolution_strict_TAGs.png")

plot_dir <- "plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ========================
# Read gene lists
# ========================
tags       <- readLines(tags_file)
non_tags   <- readLines(non_tags_file)
singletons <- readLines(singletons_file)

# ========================
# Read and filter data
# ========================
data <- read.table(input_file, header = TRUE, sep = "\t")
data <- data %>% 
  filter(is.finite(dS) & is.finite(dS_SE) & dS > 0 & dS < 6 & dS_SE < 5) %>%
  select(pair, dS)

# ========================
# Assign category based on gene IDs
# ========================
data <- data %>%
  separate(pair, into = c("gene1", "gene2"), sep = "_") %>%
  mutate(
    category = case_when(
      gene1 %in% tags & gene2 %in% tags ~ "TAG-TAG",
      (gene1 %in% tags & !(gene2 %in% tags)) | (gene2 %in% tags & !(gene1 %in% tags)) ~ "TAG-NonTAG",
      gene1 %in% non_tags & gene2 %in% non_tags ~ "NonTAG-NonTAG",
      gene1 %in% singletons | gene2 %in% singletons ~ "Singleton",
      TRUE ~ "Other"
    )
  )

# ========================
# Histogram parameters
# ========================
binwidth <- 0.01
breaks <- seq(
  from = floor(min(data$dS) / binwidth) * binwidth,
  to   = ceiling(max(data$dS) / binwidth) * binwidth,
  by   = binwidth
)

# ========================
# Plot histogram
# ========================
p <- ggplot(data, aes(x = dS, fill = category)) +
  geom_histogram(breaks = breaks, position = "stack", color = "black") +
  scale_fill_manual(
    values = c(
      "TAG-TAG" = "steelblue",
      "TAG-NonTAG" = "darkgreen",
      "NonTAG-NonTAG" = "orange",
      "Singleton" = "grey70",
      "Other" = "red"
    )
  ) +
  labs(
    title = "(H) O. brachyantha",
    x = "dS",
    y = "Number of gene pairs",
    fill = "Category"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

# ========================
# Save plot
# ========================
plot(p)
ggsave(file.path(plot_dir, output_file), plot = p, width = 7, height = 5)
message("Plot saved to: ", file.path(plot_dir, output_file))
                                     