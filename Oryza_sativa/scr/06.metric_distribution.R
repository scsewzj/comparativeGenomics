
library(Hmisc)
library(ggplot2)
library(dplyr)

# set working directory
setwd("/home/ngoc/comparativeGenomics/Oryza_sativa/")


# read filter-0 file
f0 <- read.csv("data/osativa_blastp_results_filtered.txt", sep="\t", header=F)
View(head(f0))


# setting column names
colnames(f0) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                  "query_length", "subject_length")
f0$query_coverage <- (f0$qend - f0$qstart + 1) / f0$query_length * 100
f0$subj_coverage <- (f0$send - f0$sstart + 1) / f0$subject_length * 100


max(f0$evalue)
quantiles <- quantile(f0$evalue, probs=c(0.5, 0.8, 0.95))
ggplot(f0, aes(x=evalue)) + 
  geom_histogram(bins=50) +
  geom_histogram(bins=50, fill="#d6eaf3", color="black") + 
  # geom_vline(xintercept=quantiles, linetype="dashed", color="red3") +
  # geom_text(data=data.frame(x=quantiles, y=0, label=names(quantiles)),
  #           aes(x=x, y=y, label=label), angle=0, vjust=-27, hjust=-0.2, color="red3", size=2) + 
  scale_y_continuous(expand=expansion(mult=c(0.01, 0.1))) + 
  scale_x_continuous(expand=expansion(mult=c(0.01, 0.01))) + 
  theme_minimal() +
  theme(panel.border= element_blank(), # element_rect(colour="black", fill=NA, linewidth=1),
        # axis.line.x= element_blank(), axis.line.y= element_blank(),
        plot.title=element_text(color="black", face="bold", hjust=0.5, size=13),
        panel.grid.major=element_line(color="grey97"),panel.grid.minor=element_blank(),
        axis.text=element_text(size=11, color="black"),
        axis.text.x=element_text(size=11, color="black", angle=0, hjust=1),
        legend.title=element_text(size=11, color="black", face="bold"),
        legend.text=element_text(size=11, color="black")) + 
  labs(title=" ", x="E-value", y="# Hits")



### ----- Percent identity
quantile(f0$pident)
sum(f0$pident > 30) / dim(f0)[1] * 100 # 64.12% of hits
sum(f0$pident > 50) / dim(f0)[1] * 100 # 6.97% of hits
ggplot(f0, aes(x=pident)) + 
  geom_histogram(bins=50, fill="#D2D0A0") + 
  geom_vline(xintercept=c(31, 51), linetype="dashed", color="black", size=0.6) +
  geom_text(data=data.frame(x=c(31, 51), y=0, label=c("30%", "50%")), # names(quantiles)
            aes(x=x, y=Inf, label=label), angle=0, vjust=2.0, hjust=-0.3, color="red3", size=4) +
  scale_y_continuous(expand=expansion(mult=c(0.0, 0.05))) + 
  scale_x_continuous(expand=expansion(mult=c(0.01, 0.01))) + 
  theme_minimal() +
  theme(panel.border= element_blank(), # element_rect(colour="black", fill=NA, linewidth=1),
        axis.line.x=element_line(color="black"),
        axis.line.y=element_line(color="black"),
        plot.title=element_text(color="black", face="bold", hjust=0.5, size=13),
        panel.grid.major=element_line(color="grey97"),panel.grid.minor=element_blank(),
        axis.title=element_text(size=12, color="black", face="bold"),
        axis.text=element_text(size=12, color="black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=1),
        legend.title=element_text(size=12, color="black", face="bold"),
        legend.text=element_text(size=12, color="black")) + 
  labs(title=" ", x="Percent identity (%)", y="# Hits")


### ----- Coverage
max(f0$query_coverage) # 1
max(f0$subj_coverage) # 1
quantile(f0$query_coverage, probs=c(0.3, 0.5, 0.9))
# 30%        50%        80% 
# 0.05472264 0.08588957 0.16196447 
quantile(f0$subj_coverage, probs=c(0.3, 0.5, 0.9))
# 30%        50%        80% 
# 0.05180788 0.07972603 0.15214970 

quantiles <- quantile(f0$query_coverage, probs=c(0.6, 0.8))
ggplot(f0, aes(x=query_coverage)) + 
  geom_histogram(bins=30, fill="#D2D0A0") + 
  geom_vline(xintercept=c(50, 70), linetype="dashed", color="black") +
  geom_text(data=data.frame(x=c(50, 70), y=0, label=c("50%", "70%")),
            aes(x=x, y=Inf, label=label), angle=0, vjust=2.0, hjust=-0.2, color="red3", size=4) +
  scale_y_continuous(expand=expansion(mult=c(0.0, 0.1))) + 
  scale_x_continuous(expand=expansion(mult=c(0.01, 0.01))) + 
  theme_minimal() +
  theme(panel.border= element_blank(), # element_rect(colour="black", fill=NA, linewidth=1),
        axis.line.x=element_line(color="black"),
        axis.line.y=element_line(color="black"),
        plot.title=element_text(color="black", face="bold", hjust=0.5, size=13),
        panel.grid.major=element_line(color="grey97"),panel.grid.minor=element_blank(),
        axis.title=element_text(size=12, color="black", face="bold"),
        axis.text=element_text(size=12, color="black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=1),
        legend.title=element_text(size=12, color="black", face="bold"),
        legend.text=element_text(size=12, color="black")) + 
  labs(title=" ", x="Query coverage (%)", y="# Hits")



quantiles <- quantile(f0$subj_coverage, probs=c(0.5, 0.7))
ggplot(f0, aes(x=subj_coverage)) + 
  geom_histogram(bins=30, fill="#D2D0A0") + 
  geom_vline(xintercept=c(50, 70), linetype="dashed", color="black") +
  geom_text(data=data.frame(x=c(50, 70), y=0, label=c("50%", "70%")),
            aes(x=x, y=Inf, label=label), angle=0, vjust=2.0, hjust=-0.2, color="red3", size=4) +
  scale_y_continuous(expand=expansion(mult=c(0.0, 0.1))) + 
  scale_x_continuous(expand=expansion(mult=c(0.01, 0.01))) + 
  theme_minimal() +
  theme(panel.border= element_blank(), # element_rect(colour="black", fill=NA, linewidth=1),
        axis.line.x=element_line(color="black"),
        axis.line.y=element_line(color="black"),
        plot.title=element_text(color="black", face="bold", hjust=0.5, size=13),
        panel.grid.major=element_line(color="grey97"),panel.grid.minor=element_blank(),
        axis.title=element_text(size=12, color="black", face="bold"),
        axis.text=element_text(size=12, color="black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=1),
        legend.title=element_text(size=12, color="black", face="bold"),
        legend.text=element_text(size=12, color="black")) + 
  labs(title=" ", x="Subject coverage (%)", y="# Hits")


### ----- Stats
stats <- data.frame(Sets=c(rep("L", 3), # \n(pident>0.5,\ncoverage>0.7,\ne<1e-10)
                           rep("H", 3)),
                    Stats=rep(c("hits", "seqs", "singletons"), 2),
                    Counts=c(159643, 20590, 7415, 
                             18893, 12391, 15614))
stats$Sets <- factor(stats$Sets, levels=c("L", 
                                          "H"))
stats$Stats <- factor(stats$Stats, levels=c("hits", "seqs", "singletons"))

ggplot(stats, aes(x = Sets, y = Counts)) +
  geom_col(aes(fill = Stats), position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = comma(Counts), group = Stats),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 4) +
  scale_y_continuous(limits = c(0, max(stats$Counts)*1.05),
                     expand = c(0, 0)) + # labels = comma, 
  scale_y_break(c(30000, 140000), scales = 0.5, ticklabels = c(30000, 150000)) +
  scale_fill_manual(values = c("#6caed6", "#cbaa7a", "grey")) +
  labs(title="Number of BLAST hits, sequences,\nand singletons for each dataset", x="", y="") +
  # scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 0)) +
  theme_custom() + 
  theme(
  #      panel.border=element_rect(colour="black", fill=NA, linewidth=1), #
         axis.line.x= element_line(colour="black", linewidth=0.5),
         axis.line.y= element_line(colour="black", linewidth=0.5),
         axis.text.y.right  = element_blank(),
         axis.ticks.y.right = element_blank(),
         axis.line.y.right  = element_blank(),
        axis.text.x=element_text(size=12, color="black", face="bold"),
        legend.title=element_blank(),#element_text(size=11, color="black", face="bold"),
        legend.text=element_text(size=12, color="black"),
        legend.position = "bottom")


