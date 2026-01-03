

dt.H.clean <- dt.H.clean[,c(3,2,11,12)]
dt.L.clean <- dt.L.clean[,c(3,2,11,12)]
dt.H.clean <- dt.H.clean%>%relocate(seq1, .before = seq2)
dt.L.clean <- dt.L.clean%>%relocate(seq1, .before = seq2)
dt.H.clean <- dt.H.clean[order(dt.H.clean$seq1),]
dt.L.clean <- dt.L.clean[order(dt.L.clean$seq1),]
rownames(dt.L.clean) <- NULL
rownames(dt.H.clean) <- NULL



## ----- Get the positions from the protein_info.txt

pr.info <- read.csv("data/protein_info.txt", sep="\t")
rownames(pr.info) <- pr.info$protein_id

d <- merge(dt.H.clean, pr.info[,c(5,7,8,9,12)], by.x="seq1", by.y="transcript", all.x=T)
dt.H.clean <- d
colnames(dt.H.clean)[c(5,6,7,8)] <- c("seq1_chr", "seq1_start", "seq1_end", "seq1_strand")
d <- merge(dt.H.clean, pr.info[,c(5,7,8,9,12)], by.x="seq2", by.y="transcript", all.x=T)
dt.H.clean <- d
colnames(dt.H.clean)[c(9,10,11,12)] <- c("seq2_chr", "seq2_start", "seq2_end", "seq2_strand")

d <- merge(dt.L.clean, pr.info[,c(5,7,8,9,12)], by.x="seq1", by.y="transcript", all.x=T)
dt.L.clean <- d
colnames(dt.L.clean)[c(5,6,7,8)] <- c("seq1_chr", "seq1_start", "seq1_end", "seq1_strand")
d <- merge(dt.L.clean, pr.info[,c(5,7,8,9,12)], by.x="seq2", by.y="transcript", all.x=T)
dt.L.clean <- d
colnames(dt.L.clean)[c(9,10,11,12)] <- c("seq2_chr", "seq2_start", "seq2_end", "seq2_strand")
rm(d)



## ----- Stratifying
d <- unique(unlist(dt.H.clean%>%filter(dS<0.18)%>%select(c(seq1, seq2))))
writeLines(d, "paml-workflow/pantherdb/high/input/dS.lt0.18.txt") # 1657g, 1518 pairs

d <- unique(unlist(dt.H.clean%>%filter(dS<1.8 & dS>1.3)%>%select(c(seq1, seq2))))
writeLines(d, "paml-workflow/pantherdb/high/input/dS.gt1.3_lt1.8.txt") # 3218g, 3741 pairs

d <- unique(unlist(dt.H.clean%>%filter(dS>3.5)%>%select(c(seq1, seq2))))
writeLines(d, "paml-workflow/pantherdb/high/input/dS.gt3.5.txt") # 3460g, 7547 pairs

rm(d)



## ----- Plotting

list.files("paml-workflow/pantherdb/high/output/")
go.df <- load_go_results("paml-workflow/pantherdb/high/output/Slim-dS.gt3.5.txt", fdr_filter=0.05)
# width=573 & height=461
ggplot(go.df[1:15,], aes(x=fold_enrichment, y=reorder(GO_term, fold_enrichment, decreasing=F))) + 
  geom_point(aes(fill=FDR, size=number_in_list), shape=21, color="black") +
  scale_y_discrete(expand=expansion(mult=c(0.05, 0.05)), labels=label_wrap(53)) +
  scale_x_continuous(expand=expansion(mult=c(0.2, 0.2))) + 
  scale_fill_continuous(high="#DAE2DD", low="#4A6F56", name="FDR", 
                        limits=c(0, 0.05), breaks=c(0.01, 0.02, 0.03, 0.04)) +
  scale_size(range=c(2, 8), name="Number in list", 
             breaks=c(5, 10, 20, 40, 80), limits=c(2, 81)) +
  labs(x="Fold enrichment", y=NULL, 
       title="PantherDB SLIM-BP\ndS>3.5 (3460 genes, 7547 pairs)") +
  guides(size=guide_legend(order=1, theme=theme(legend.text=element_text(size=10, color="black"),
                                                legend.title=element_text(size=10, face="bold", color="black"))),
         color=guide_colorbar(order=2, theme=theme(legend.text=element_text(size=10, color="black"),
                                                   legend.title=element_text(size=10, face="bold", color="black")))) + 
  theme_minimal() +
  theme_custom()



# ---- clean a bit
ds.H <- dt.H.clean %>%
  transmute(g1 = pmin(seq1, seq2),
            g2 = pmax(seq1, seq2),
            dS = dS)%>%
  distinct()

ds.L <- dt.L.clean %>%
  transmute(g1 = pmin(seq1, seq2),
            g2 = pmax(seq1, seq2),
            dS = dS)%>%
  distinct()

























