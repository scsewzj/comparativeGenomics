library(ggplot2)

# Read the numeric .gene_type file
gene_type <- read.table("./ob.gene_type", header = FALSE,
                        col.names = c("Gene", "TypeCode"))

# Map numeric codes to descriptive labels
dup_map <- c("0"="Singleton",
             "4"="WGD",
             "3"="Tandem",
             "2"="Proximal",
             "1"="Dispersed")

gene_type$Type <- dup_map[as.character(gene_type$TypeCode)]

# View
head(gene_type)

# Read MCScanX duplication result
dup <- gene_type
# head(dup)

# Count each type
dup_count <- as.data.frame(table(dup$Type))
colnames(dup_count) <- c("Type", "Count")
# head(dup_count)

# Plot
ggplot(dup_count, aes(x = Type, y = Count)) +
  geom_bar(stat = "identity",  fill = "lightgreen", color = "black") +
  theme_minimal() +
  labs(
    title = "Distribution of Gene Duplication Types",
    x = "Duplication Type",
    y = "Number of Genes"
  )
