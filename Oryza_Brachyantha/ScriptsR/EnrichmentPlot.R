# Load libraries
library(ggplot2)
library(dplyr)
library(readr)


read_panther_ora <- function(file_path) {
  # 1. Read first 50 lines to detect header
  lines <- readLines(file_path, n = 50)
  
  # Find the line number containing "Fold Enrichment" (header line)
  header_line_num <- grep("fold Enrichment", lines, ignore.case = TRUE)
  if(length(header_line_num) == 0) stop("Cannot find header line containing 'fold Enrichment'")
  
  # 2. Read the table starting from the header line
  data <- read_tsv(file_path, skip = header_line_num - 1, col_names = TRUE, na = c("", "NA"))
  
  # 3. Identify important columns (by approximate name matching)
  colnames(data) <- gsub("\\s+", " ", colnames(data)) # clean extra spaces
  go_col <- grep("GO", colnames(data), ignore.case = TRUE)[1]
  fold_col <- grep("fold Enrichment", colnames(data), ignore.case = TRUE)[1]
  fdr_col <- grep("FDR", colnames(data), ignore.case = TRUE)[1]
  
  if(any(is.na(c(go_col, fold_col, fdr_col)))) stop("Cannot detect required columns")
  
  # 4. Keep only these columns
  data_clean <- data %>%
    select(GO_term = go_col,
           FoldEnrichment = fold_col,
           FDR = fdr_col) %>%
    filter(!is.na(FoldEnrichment))  # remove empty rows
  
  # 5. Convert to numeric
  data_clean <- data_clean %>%
    mutate(FoldEnrichment = as.numeric(FoldEnrichment),
           FDR = as.numeric(FDR))
  
  return(data_clean)
}

# Usage example

# Read your PANTHER file (tab-separated)
file_path <- "../PANTHER/H-analysis-OverR-MolF.txt"
df <- read_panther_ora(file_path)
head(df)

PLOT_DIR <- "plots"
OUTPUT_FILE <- "H-EnrichmentFold.png"

# Plot
p <- ggplot(df, aes(x = FoldEnrichment, y = GO_term)) +
  geom_point(aes(size = -log10(FDR), color = -log10(FDR))) +
  scale_color_gradient(low = "skyblue", high = "red") +
  labs(x = "Fold Enrichment", y = "", 
       color = "-log10(FDR)", size = "-log10(FDR)",
       title = "GO Slim Biological Process Enrichment") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))


plot(p)

ggsave(file.path(PLOT_DIR, OUTPUT_FILE), plot = p, width = 7, height = 5)
message("Plot saved to: ", file.path(PLOT_DIR, OUTPUT_FILE))