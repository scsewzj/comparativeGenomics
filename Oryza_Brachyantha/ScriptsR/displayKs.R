
# Load libraries
library(ggplot2)

## Simple R Script: Histogram + Density
# Read your Ks values (one column file)
# Assumes file has a column named "Ks"
ks <- read.table("Ks_values.txt", header=TRUE)

# Basic histogram + density
ggplot(ks, aes(x = Ks)) +
  geom_histogram(aes(y = ..density..), bins = 50) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(title = "Ks Distribution",
       x = "Ks",
       y = "Density")
       
## Density-only plot (classic Ks peak plot)
ggplot(ks, aes(x = Ks)) +
  geom_density(fill = NA) +
  theme_minimal() +
  labs(title = "Ks Density Plot",
       x = "Ks",
       y = "Density")
       
## Multiple clusters
ggplot(ks, aes(x = Ks)) +
  geom_density() +
  facet_wrap(~ Cluster, scales = "free_y") +
  theme_minimal() +
  labs(title = "Ks Density by Cluster",
       x = "Ks",
       y = "Density")
       
## Combined violin + boxplot
ggplot(ks, aes(x = Cluster, y = Ks)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme_minimal() +
  labs(title = "Ks Distribution Across Clusters")
  

# ---------------------------------------------
# ------------ ks_analysis.R ------------------
# ---------------------------------------------
library(ggplot2)
library(dplyr)

# Load Ks data
ks <- read.table("Ks_all_clusters.tsv", header = TRUE, sep = "\t")

# ---- Global density plot ----
p1 <- ggplot(ks, aes(Ks)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(title = "Global Ks Density", x = "Ks", y = "Density")

ggsave("Ks_density_global.pdf", p1, width = 7, height = 5)

# ---- Peak detection ----
dens <- density(ks$Ks, na.rm = TRUE)
x <- dens$x
y <- dens$y

# Find local maxima
peaks <- which(diff(sign(diff(y))) == -2) + 1
peak_positions <- x[peaks]
peak_values <- y[peaks]

# Save peaks
peak_df <- data.frame(Ks_peak = peak_positions, Density = peak_values)
write.table(peak_df, "Ks_detected_peaks.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Plot with peaks labeled ----
p2 <- ggplot() +
  geom_line(aes(x, y), color = 1) +
  geom_point(aes(peak_positions, peak_values), color = 2, size = 3) +
  geom_text(aes(peak_positions, peak_values, label = round(peak_positions, 3)),
            vjust = -1) +
  theme_minimal() +
  labs(title="Ks Density With Detected Peaks", x="Ks", y="Density")

ggsave("Ks_density_peaks.pdf", p2, width = 7, height = 5)

# -----------Pairwise comparison plots ------------ #
# Density by cluster (faceted Ks)
p3 <- ggplot(ks, aes(Ks)) +
  geom_density() +
  facet_wrap(~ Cluster, scales = "free_y") +
  theme_minimal() +
  labs(title="Ks Density by Cluster")
ggsave("Ks_density_by_cluster.pdf", p3, width = 10, height = 8)

# Violin + boxplot across clusters
p4 <- ggplot(ks, aes(x = Cluster, y = Ks)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  theme_minimal() +
  labs(title="Ks by Cluster (Violin + Boxplot)")
ggsave("Ks_violin_by_cluster.pdf", p4, width = 12, height = 6)
