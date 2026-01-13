# Load libraries
library(ggplot2)
library(dplyr)
library(readr)


read_panther_ora <- function(file_path) {
  library(readr)
  library(dplyr)
  
  # 1. Read first 50 lines to detect header
  lines <- readLines(file_path, n = 50)
  
  # Find the line number containing "Fold Enrichment" (header line)
  header_line_num <- grep("fold Enrichment", lines, ignore.case = TRUE)
  if(length(header_line_num) == 0) stop("Cannot find header line containing 'fold Enrichment'")
  
  # 2. Read the table starting from the header line
  data <- read_tsv(file_path, skip = header_line_num - 1, col_names = TRUE, na = c("", "NA"))
  
  # 3. Clean column names
  colnames(data) <- gsub("\\s+", " ", colnames(data)) # clean extra spaces
  
  # 4. Identify important columns
  go_col <- grep("GO", colnames(data), ignore.case = TRUE)[1]
  genes_col <- 3  # the 3rd column in your table, number of genes in your list
  fold_col <- grep("fold Enrichment", colnames(data), ignore.case = TRUE)[1]
  fdr_col <- grep("FDR", colnames(data), ignore.case = TRUE)[1]
  
  if(any(is.na(c(go_col, fold_col, fdr_col)))) stop("Cannot detect required columns")
  
  # 5. Keep only the relevant columns
  data_clean <- data %>%
    select(
      GO_term = go_col,
      GeneCount = genes_col,        # new column
      FoldEnrichment = fold_col,
      FDR = fdr_col
    ) %>%
    filter(!is.na(FoldEnrichment))  # remove empty rows
  
  # 6. Convert numeric columns
  data_clean <- data_clean %>%
    mutate(
      GeneCount = as.numeric(GeneCount),
      FoldEnrichment = as.numeric(FoldEnrichment),
      FDR = as.numeric(FDR)
    )
  
  # 7. Clean and filter
  data_clean <- data_clean %>%
    mutate(GO_term = gsub("\\s*\\(GO:[0-9]+\\)", "", GO_term)) %>% 
    mutate(GO_term = trimws(GO_term)) %>%
    filter(
      FDR < 0.05,
      !grepl("unclassified", GO_term, ignore.case = TRUE),
      !grepl("biological_process", GO_term, ignore.case = TRUE),
      !grepl("molecular_function", GO_term, ignore.case = TRUE)
    ) %>%
    arrange(FDR) %>%      # sort by lowest FDR
    slice_head(n = 15)    # keep top 15
  
  return(data_clean)
}


# Usage example

# Read your PANTHER file (tab-separated)
file_path_singleton <- "../PANTHER/H_Analysis/H-single-analysis-MF.txt"
file_path_nonTag <- "../PANTHER/H_Analysis/H-nonTAGs-analysis-MF.txt"
file_path_Tag <- "../PANTHER/H_Analysis/H-TAGs-analysis-MF.txt"

# ------------------------
# 1. Singletons
# ------------------------
df <- read_panther_ora(file_path_singleton)
head(df)

PLOT_DIR <- "../plots"
OUTPUT_FILE <- "H-single-EnrichmentFold-MolF.png"

# ------------------------
# Create output directory
# ------------------------
if (!dir.exists(PLOT_DIR)) {
  dir.create(PLOT_DIR, recursive = TRUE)
}

# Plot
p <- ggplot(df, aes(x = FoldEnrichment, y = GO_term)) +
  geom_point(aes(size = GeneCount, color =  FDR)) +  # size = GeneCount
  scale_color_gradient(low = "skyblue", high = "red") +
  labs(
    x = "Fold Enrichment", 
    y = "", 
    color = " FDR", 
    size = "Genes in list",  # updated legend label
    title = "(H) Singletons"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

plot(p)

ggsave(file.path(PLOT_DIR, OUTPUT_FILE), plot = p)
message("Plot saved to: ", file.path(PLOT_DIR, OUTPUT_FILE))

# ------------------------
# 2. non-TAG
# ------------------------
df <- read_panther_ora(file_path_nonTag)
OUTPUT_FILE <- "H-nonTAGs-EnrichmentFold-MolF.png"

# Plot
p <- ggplot(df, aes(x = FoldEnrichment, y = GO_term)) +
  geom_point(aes(size = GeneCount, color =  FDR)) +
  scale_color_gradient(low = "skyblue", high = "red") +
  labs(x = "Fold Enrichment", y = "", 
       color = "FDR", size = "Genes in list",
       title = "(H) non-TAGs") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))


plot(p)

ggsave(file.path(PLOT_DIR, OUTPUT_FILE), plot = p)
message("Plot saved to: ", file.path(PLOT_DIR, OUTPUT_FILE))


# ------------------------
# 3. TAGs
# ------------------------
df <- read_panther_ora(file_path_Tag)
OUTPUT_FILE <- "H-TAGs-EnrichmentFold-MolF.png"

# Plot
p <- ggplot(df, aes(x = FoldEnrichment, y = GO_term)) +
  geom_point(aes(size = GeneCount, color =  FDR)) +
  scale_color_gradient(low = "skyblue", high = "red") +
  labs(x = "Fold Enrichment", y = "", 
       color = "FDR", size = "Genes in list",
       title = "(H) TAGs: 0 spacers") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))


plot(p)

ggsave(file.path(PLOT_DIR, OUTPUT_FILE), plot = p)
message("Plot saved to: ", file.path(PLOT_DIR, OUTPUT_FILE))