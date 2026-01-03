
library(ggbreak)
library(grid)
library(scales)

setwd("/home/ngoc/comparativeGenomics/Oryza_sativa/")


theme_custom <- function() {
  theme_minimal() +
    theme(panel.border= element_blank(), # element_rect(colour="black", fill=NA, linewidth=1),
          # axis.line.x= element_blank(), axis.line.y= element_blank(),
          plot.title=element_text(color="black", face="bold", hjust=0.5, size=13),
          panel.grid.major=element_line(color="grey97"),panel.grid.minor=element_blank(),
          axis.text=element_text(size=11, color="black"),
          axis.text.x=element_text(size=11, color="black", angle=0, hjust=1),
          axis.text.y.right  = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.title.y.right = element_blank(),
          legend.title=element_text(size=11, color="black", face="bold"),
          legend.text=element_text(size=11, color="black")) 
}

high <- read.csv("mcl/output/H.tabular", sep="\t", header=F)
low <- read.csv("mcl/output/L.tabular", sep="\t", header=F)

# calculate the size of each cluster
stats.high <- sapply(1:nrow(high), function(i) sum(high[i, ] != "" & !is.na(high[i, ])))
stats.low <- sapply(1:nrow(low), function(i) sum(low[i, ] != "" & !is.na(low[i, ])))

summary(stats.high)
summary(stats.low)

# size of genes
high.genes <- as.vector(as.matrix(high)) 
high.genes <- unique(high.genes[!is.na(high.genes) & high.genes != ""])
low.genes <- as.vector(as.matrix(low)) 
low.genes <- unique(low.genes[!is.na(low.genes) & low.genes != ""])
length(high.genes) # 12391
length(low.genes) # 20590 # all input sequences are present here.


# distribution of cluster size
ggplot(as.data.frame(stats.high), aes(x=stats.high)) + 
  geom_histogram(bins=30, fill="#bf828c", color="black") + 
  scale_y_break(c(1000, 2500), scales = 0.5, ticklabels = c(1000, 2600)) +
  scale_y_continuous(expand=expansion(mult=c(0.01, 0.1))) + 
  scale_x_continuous(expand=expansion(mult=c(0.05, 0.01))) + 
  labs(title="High", x="Cluster size", y="Frequency") + theme_custom() 


ggplot(as.data.frame(stats.low), aes(x=stats.low)) + 
  geom_histogram(bins=30, fill="#9f9ac7", color="black") + 
  scale_y_continuous(expand=expansion(mult=c(0.01, 0.1))) + 
  scale_x_continuous(expand=expansion(mult=c(0.05, 0.01))) + 
  scale_y_break(c(750, 3500), scales = 0.5, ticklabels = c(750, 3800)) +
  labs(title="Low", x="Cluster size", y="Frequency") + 
  theme_custom()



## ----- Statistics
stats.high.df <- data.frame(cluster_index=seq(1, length(stats.high), 1), size = stats.high)
stats.low.df <- data.frame(cluster_index=seq(1, length(stats.low), 1), size = stats.low)

stats.summary <- data.frame(
  Cluster_size = c("all", "max size", "min size", "mean+-SD size", "2", "3–9", "10–19", ">=20"),
  High = c(length(stats.high),
           max(stats.high),
           min(stats.high),
           paste0(mean(stats.high), "+-", sd(stats.high)),
           sum(stats.high == 2),
           sum(stats.high >= 3 & stats.high <= 9),
           sum(stats.high >= 10 & stats.high <= 19),
           # sum(stats.high >= 3 & stats.high < 6),
           # sum(stats.high >= 6 & stats.high < 11),
           # sum(stats.high >= 11 & stats.high < 16),
           # sum(stats.high >= 16 & stats.high < 20),
           sum(stats.high >= 20)),
  Low = c(length(stats.low),
          max(stats.low),
          min(stats.low),
          paste0(mean(stats.low), "+-", sd(stats.low)),
          sum(stats.low == 2),
          sum(stats.low >= 3 & stats.low <= 9),
          sum(stats.low >= 10 & stats.low <= 19),
          # sum(stats.low >= 3 & stats.low < 6),
          # sum(stats.low >= 6 & stats.low < 11),
          # sum(stats.low >= 11 & stats.low < 16),
          # sum(stats.low >= 16 & stats.low < 20),
          sum(stats.low >= 20)),
  row.names = NULL)


## ----- extracting clusters of different sizes for H
high.genes.gt20 <- as.vector(unlist(high[which(stats.high>=20),]))
high.genes.gt20 <- unique(high.genes.gt20[!is.na(high.genes.gt20) & high.genes.gt20 != ""])

high.genes.gt10.lt20 <- as.vector(unlist(high[which(stats.high<20 & stats.high>=10),]))
high.genes.gt10.lt20 <- unique(high.genes.gt10.lt20[!is.na(high.genes.gt10.lt20) & high.genes.gt10.lt20 != ""])

high.genes.gt3.lt10 <- as.vector(unlist(high[which(stats.high<10 & stats.high>=3),]))
high.genes.gt3.lt10 <- unique(high.genes.gt3.lt10[!is.na(high.genes.gt3.lt10) & high.genes.gt3.lt10 != ""])

high.genes.eq2 <- as.vector(unlist(high[which(stats.high==2),]))
high.genes.eq2 <- unique(high.genes.eq2[!is.na(high.genes.eq2) & high.genes.eq2 != ""])


## ----- extracting clusters of different sizes for L
low.genes.gt20 <- as.vector(unlist(low[which(stats.low>=20),]))
low.genes.gt20 <- unique(low.genes.gt20[!is.na(low.genes.gt20) & low.genes.gt20 != ""])

low.genes.gt10.lt20 <- as.vector(unlist(low[which(stats.low<20 & stats.low>=10),]))
low.genes.gt10.lt20 <- unique(low.genes.gt10.lt20[!is.na(low.genes.gt10.lt20) & low.genes.gt10.lt20 != ""])

low.genes.gt3.lt10 <- as.vector(unlist(low[which(stats.low<10 & stats.low>=3),]))
low.genes.gt3.lt10 <- unique(low.genes.gt3.lt10[!is.na(low.genes.gt3.lt10) & low.genes.gt3.lt10 != ""])

low.genes.eq2 <- as.vector(unlist(low[which(stats.low==2),]))
low.genes.eq2 <- unique(low.genes.eq2[!is.na(low.genes.eq2) & low.genes.eq2 != ""])


writeLines(low.genes.gt20, "mcl/pantherdb/input/low.genes.gt20.txt")
writeLines(low.genes.gt10.lt20, "mcl/pantherdb/input/low.genes.gt10.lt20.txt")
writeLines(low.genes.gt3.lt10, "mcl/pantherdb/input/low.genes.gt3.lt10.txt")
writeLines(low.genes.eq2, "mcl/pantherdb/input/low.genes.eq2.txt")
writeLines(high.genes.gt20, "mcl/pantherdb/input/high.genes.gt20.txt")
writeLines(high.genes.gt10.lt20, "mcl/pantherdb/input/high.genes.gt10.lt20.txt")
writeLines(high.genes.gt3.lt10, "mcl/pantherdb/input/high.genes.gt3.lt10.txt")
writeLines(high.genes.eq2, "mcl/pantherdb/input/high.genes.eq2.txt")


## ----- Plot
load_go_results <- function(path, plus_only=T, fdr_filter=0.05) {
  df <- read.delim(path, header = TRUE, stringsAsFactors = FALSE)
  colnames(df) <- c("GO", "number_in_ref", "number_in_list", "expected", "plus_minus", "fold_enrichment", "pval", "FDR")
  df$fold_enrichment <- as.numeric(df$fold_enrichment)
  df <- df%>%filter(!is.na(fold_enrichment))
  df$GO_term <- sub("\\s*\\(GO:\\d+\\)$", "", df$GO)
  df$GOID <- sub(".*\\((GO:\\d+)\\)$", "\\1", df$GO)
  df <- dplyr::relocate(df, GOID, GO_term, .before = number_in_ref)
  if (plus_only) {
    df <- dplyr::filter(df, plus_minus == "+")
  }
  df <- dplyr::filter(df, FDR<fdr_filter)
  df <- df%>%filter(GOID!="Unclassified (UNCLASSIFIED)")
  df[,-1]
}

list.files("mcl/pantherdb/output/")
go.df <- load_go_results("mcl/pantherdb/output/Slim-low_TAGs_0.txt", fdr_filter=0.05)
# width=573 & height=461
ggplot(go.df[c(1:15),], aes(x=fold_enrichment, y=reorder(GO_term, fold_enrichment, decreasing=F))) + 
  geom_point(aes(fill=FDR, size=number_in_list), shape=21, color="black") +
  scale_y_discrete(expand=expansion(mult=c(0.05, 0.05)), labels=label_wrap(60)) +
  scale_x_continuous(expand=expansion(mult=c(0.2, 0.2))) + 
  scale_fill_continuous(high="#DAE2DD", low="#4A6F56", name="FDR", 
                        limits=c(0, 0.05), breaks=c(0.01, 0.02, 0.03, 0.04)) +
  scale_size(range=c(2, 8), name="Number in list", 
             breaks=c(5, 10, 20, 40, 80), limits=c(1, 80)) +
  labs(x="Fold enrichment", y=NULL, 
       title="(L) TAGs:0 spacers") +
  guides(size=guide_legend(order=1, theme=theme(legend.text=element_text(size=10, color="black"),
                                                legend.title=element_text(size=10, face="bold", color="black"))),
         color=guide_colorbar(order=2, theme=theme(legend.text=element_text(size=10, color="black"),
                                                   legend.title=element_text(size=10, face="bold", color="black")))) + 
  theme_minimal() +
  theme_custom()






## ----- Helper functions
theme_custom <- function() {
  theme_minimal() +
    theme(panel.border=element_rect(color="grey80", fill=NA, linewidth=1),
          axis.line.x=element_blank(),
          axis.line.y=element_blank(),
          plot.title=element_text(size=10, color="black", face="bold", hjust=0.5),
          panel.grid.major=element_line(color="grey97"),
          panel.grid.minor=element_blank(),
          axis.text=element_text(size=10, color="black"),
          axis.text.x=element_text(size=10, color="black", hjust=1),
          axis.title.x=element_text(size=10, color="black", face="bold", margin=margin(t=10)),
          axis.title.y=element_blank(),
          legend.text=element_text(size=10, color="black"),
          legend.title=element_text(size=10, face="bold", color="black")) 
}















