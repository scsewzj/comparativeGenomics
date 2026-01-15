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

L_file  <- ifelse(length(args) >= 1, args[1], "../KaKs_all_ks_input_moderate.tsv")
H_file  <- ifelse(length(args) >= 1, args[1], "../KaKs_all_ks_input_strict.tsv")
output_file <- ifelse(length(args) >= 2, args[2], "dS_SE_vs_dS.png")

plot_dir <- "plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ========================
# Load data
# ========================
missing_colsL <- setdiff(required_cols, colnames(dfL))
missing_colsH <- setdiff(required_cols, colnames(dfH))

if (length(missing_colsL) > 0) {
  stop(paste("dfL is missing required columns:",
             paste(missing_colsL, collapse = ", ")))
}

if (length(missing_colsH) > 0) {
  stop(paste("dfH is missing required columns:",
             paste(missing_colsH, collapse = ", ")))
}
# Remove invalid values
threshold_SE <- 5
threshold <- 6
df1 <- dfL[is.finite(dfL$dS) & is.finite(dfL$dS_SE) & dfL$dS >= 0 & dfL$dS < threshold,]
df2 <- dfL[is.finite(dfL$dS) & is.finite(dfL$dS_SE) & dfL$dS >= 0 & dfL$dS < threshold & dfL$dS_SE <= threshold_SE, ]
df3 <- dfH[is.finite(dfH$dS) & is.finite(dfH$dS_SE) & dfH$dS >= 0 & dfH$dS < threshold,]
df4 <- dfH[is.finite(dfH$dS) & is.finite(dfH$dS_SE) & dfH$dS >= 0 & dfH$dS < threshold & dfH$dS_SE <= threshold_SE, ]

# ========================
# Plot
# ========================
library(patchwork)

# Create individual plots
p1 <- ggplot(df1, aes(x = dS, y = dS_SE)) +
  geom_point(alpha = 0.5, size = 1) +
  labs(title = paste("A. (L) dS <6\n",dim(df1)[1]*dim(df1)[2],"pairs"),
       x = "dS", y = "SE(dS)") +
  theme_minimal(base_size = 12)

p2 <- ggplot(df2, aes(x = dS, y = dS_SE)) +
  geom_point(alpha = 0.5, size = 1) +
  labs(title = paste("B. (L) dS <6, SE< 5\n",dim(df2)[1]*dim(df2)[2],"pairs"),
       x = "dS", y = "SE(dS)") +
  theme_minimal(base_size = 12)

p3 <- ggplot(df3, aes(x = dS, y = dS_SE)) +
  geom_point(alpha = 0.5, size = 1) +
  labs(title = paste("C. (H) dS <6\n",dim(df3)[1]*dim(df3)[2],"pairs"),
       x = "dS", y = "SE(dS)") +
  theme_minimal(base_size = 12)

p4 <- ggplot(df4, aes(x = dS, y = dS_SE)) +
  geom_point(alpha = 0.5, size = 1) +
  labs(title = paste("D. (H) dS <6, SE< 5\n",dim(df4)[1]*dim(df4)[2],"pairs"), x = "dS", y = "SE(dS)") +
  theme_minimal(base_size = 12)

# Combine in 2x2 layout
combined_plot <- (p1 | p2) / (p3 | p4)

# Display
plot(combined_plot)
# ========================
# Save
# ========================
# ggsave(file.path(plot_dir, paste0("Moderate_cutoff_", threshold, "_SE_", threshold_SE, "_", output_file)), plot = p, width = 7, height = 5)

message("Plot saved to: ", file.path(plot_dir, paste0("cutoff_", threshold, "_SE_", threshold_SE, "_", output_file)))
