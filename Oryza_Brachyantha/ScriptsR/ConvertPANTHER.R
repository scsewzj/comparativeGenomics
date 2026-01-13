library(dplyr)
library(tidyr)
library(stringr)

# Load and inspect
panther <- read.table("TAGs-pantherGeneList.txt",
                      sep = "\t",
                      header = FALSE,
                      quote = "",
                      fill = TRUE)

# Keep only true GO terms
panther_go <- panther %>%
  filter(str_detect(V1, "^GO:"))

# Extract EnsemblGenome gene IDs
panther_go <- panther_go %>%
  mutate(GeneID = str_extract(V2, "Os[0-9]{2}g[0-9]+"))

# Handle rows with multiple genes
panther_go <- panther_go %>%
  separate_rows(GeneID, sep = ",")

# Final TERM2GENE table
term2gene <- panther_go %>%
  select(GO = V1, GeneID) %>%
  distinct()

GENE_FILE <- "../TAGS/TAGs_OB2GO_Overrep.txt"
your_gene_list <- read.table(file = GENE_FILE, header = F)


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")
library(ggplot2)
library(clusterProfiler)
# Run GO enrichment (FDR-controlled)
ego <- enricher(
  gene = your_gene_list$V1,     # duplicated genes
  TERM2GENE = term2gene,
  pAdjustMethod = "fdr",
  qvalueCutoff = 0.05,
  minGSSize = 1,    # allow small GO terms
  maxGSSize = 75  # allow large GO terms
)

# Dotplot: FDR = color, genes-in-list = size
dotplot(
  ego,
  showCategory = 20,
  color = "p.adjust",
  size = "Count"
) +
  scale_color_continuous(low = "red", high = "blue") +
  theme_minimal()
