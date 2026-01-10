#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
})

# ========================
# Arguments / files
# ========================
args <- commandArgs(trailingOnly = TRUE)

KaKs_strict_file   <- ifelse(length(args) >= 1, args[1], "./KaKs_all_ks_input_strict.tsv")
KaKs_moderate_file <- ifelse(length(args) >= 2, args[2], "./KaKs_all_ks_input_moderate.tsv")

moderate_dup_file  <- ifelse(length(args) >= 3, args[3], "moderate.duplicate.tsv")
moderate_singleton_file <- ifelse(length(args) >= 4, args[4], "moderate.singleton.tsv")

strict_dup_file    <- ifelse(length(args) >= 5, args[5], "strict.duplicate.tsv")
strict_singleton_file <- ifelse(length(args) >= 6, args[6], "strict.singleton.tsv")

plot_dir <- "plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ========================
# Helper function: read KaKs table
# ========================
read_kaks <- function(file) {
  df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  stopifnot(all(c("pair",	"dS",	"dN",	"omega",	"dN_SE",	"dS_SE") %in% colnames(df)))
  return(df)
}

# ========================
# Helper function: read duplicated / singleton gene list
# ========================
read_gene_list <- function(file) {
  readLines(file)
}

# ========================
# Load Ka/Ks data
# ========================
kaks_strict   <- read_kaks(KaKs_strict_file)
kaks_moderate <- read_kaks(KaKs_moderate_file)
kaks_strict <- kaks_strict[kaks_strict$dS < 6, ]
kaks_moderate <- kaks_moderate[kaks_moderate$dS < 6, ]
# Load duplicates
dup_strict     <- read_gene_list(strict_dup_file)
dup_moderate   <- read_gene_list(moderate_dup_file)

# Load singletons if needed
singleton_strict   <- read_gene_list(strict_singleton_file)
singleton_moderate <- read_gene_list(moderate_singleton_file)

# ========================
# Function: get group gene pairs
# ========================
get_gp_pairs <- function(kaks_df, gp_list) {

  kaks_df %>%
    # split "gene1_gene2" into two columns
    tidyr::separate(
      col = pair,
      into = c("gene1", "gene2"),
      sep = "_",
      remove = FALSE
    ) %>%
    # keep pairs where at least one gene is duplicated
    dplyr::filter(gene1 %in% gp_list | gene2 %in% gp_list)
}
# ========================
# Compute average omega for each stringency
# ========================
avg_omega_strict_dup <- get_gp_pairs(kaks_strict, dup_strict) %>%
  summarise(avg_omega = mean(omega, na.rm = TRUE))

avg_omega_strict_sing <- get_gp_pairs(kaks_strict, singleton_strict) %>%
  summarise(avg_omega = mean(omega, na.rm = TRUE))

avg_omega_moderate_dup <- get_gp_pairs(kaks_moderate, dup_moderate) %>%
  summarise(avg_omega = mean(omega, na.rm = TRUE))

avg_omega_moderate_sing <- get_gp_pairs(kaks_moderate, singleton_moderate) %>%
  summarise(avg_omega = mean(omega, na.rm = TRUE))

# Combine into one dataframe
df_plot <- data.frame(
  Stringency = c("Moderate", "Moderate", "Strict", "Strict"),
  Category   = c("Duplicated", "Singleton", "Duplicated", "Singleton"),
  AvgOmega   = c(
    avg_omega_moderate_dup$avg_omega,
    avg_omega_moderate_sing$avg_omega,
    avg_omega_strict_dup$avg_omega,
    avg_omega_strict_sing$avg_omega
  )
)


# ========================
# Bar chart
# ========================
p <- ggplot(df_plot, aes(x = Stringency, y = AvgOmega, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    aes(label = round(AvgOmega, 3)),
    position = position_dodge(width = 0.7),
    vjust = -0.5,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "Duplicated" = "#1b9e77",
      "Singleton"  = "#7570b3"
    )
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = expression("Average " * omega * " by Stringency and Gene Category"),
    x = "Stringency",
    y = expression("Average " * omega),
    fill = "Gene category"
  )


# ========================
# Save plot
# ========================
ggsave(file.path(plot_dir, "avg_omega_duplicates.png"), plot = p, width = 6, height = 4)
