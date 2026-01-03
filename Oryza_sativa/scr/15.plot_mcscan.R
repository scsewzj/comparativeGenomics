

tandem.dt <- read.csv("mcscan-run/osativa_H.tandem", header=F)
dim(tandem.dt) # 4225
sum(dt.H.clean$seq1_type == dt.H.clean$seq2_type & dt.H.clean$seq1_type=="Tandem") # 6108



### ---- dup types 
dup.types <- read.csv("mcscan-run/osativa_H.gene_type", sep="\t", header=F)
colnames(dup.types) <- c("id", "code")
dup.types <- dup.types%>%mutate(type=case_when(code==0~"Singleton",
                                               code==1~"Dispersed",
                                               code==2~"Proximal",
                                               code==3~"Tandem",
                                               code==4~"WGD / segmental",
                                               TRUE~"Unknown"))

table(dup.types$type)
# Dispersed       Proximal       Singleton        Tandem WGD / segmental 
# 5898            1960           24170            7290            2902 


all(unique(c(dt.H.clean$seq1, dt.H.clean$seq2)) %in% dup.types$id) # TRUE
dim(dt.H.clean) # 34920    13
dim(dup.types) # 42220     3

dt.H.clean$seq1_type <- dup.types$type[match(dt.H.clean$seq1, dup.types$id)]
dt.H.clean$seq2_type <- dup.types$type[match(dt.H.clean$seq2, dup.types$id)]

table(dt.H.clean$seq1_type)
# Dispersed        Proximal          Tandem WGD / segmental 
# 11548            4806           12262            6304 
table(dt.H.clean$seq2_type)
# Dispersed        Proximal          Tandem WGD / segmental 
# 11511            4778           12389            6242 

sum(dt.H.clean$seq1_type == dt.H.clean$seq2_type) # 18256
sum(dt.H.clean$seq1_type == dt.H.clean$seq2_type) # 16664

dt.L.clean$seq1_type <- dup.types$type[match(dt.L.clean$seq1, dup.types$id)]
dt.L.clean$seq2_type <- dup.types$type[match(dt.L.clean$seq2, dup.types$id)]

# attempting the priority:
# WGD/segmental > tandem > proximal > dispersed)
priority <- c("WGD / segmental"=1, "Tandem"=2, "Proximal"=3, "Dispersed"=4)
dt.H.clean <- dt.H.clean%>%mutate(rank1=priority[seq1_type],
                                  rank2=priority[seq2_type])
dt.H.clean <- dt.H.clean%>%mutate(pair_type=case_when(
  rank1<rank2~seq1_type, rank2<rank1~seq2_type, TRUE~seq1_type))
dt.H.clean$pair_type <- factor(dt.H.clean$pair_type,
                               levels = c("WGD / segmental","Tandem","Proximal","Dispersed"))
table(dt.H.clean$pair_type)

dup.colors <- c("WGD / segmental"="#cb7a87", "Tandem"="#e7ba51",
                "Proximal"="#78b43e", "Dispersed"="#999999")
ggplot(dt.H.clean, aes(x=dS, fill=pair_type)) + 
  geom_histogram(bins=100, alpha=0.7) + # , color="black"
  scale_fill_manual(values=dup.colors) +
  scale_y_continuous(expand=expansion(mult=c(0.0, 0.05))) + 
  scale_x_continuous(expand=expansion(mult=c(0.0, 0.01))) + 
  labs(title=bquote(bolditalic("Oryza sativa")~bold("(H)")), x="dS", y="# duplicated genes") + 
  theme_custom()


dt.L.clean <- dt.L.clean%>%mutate(rank1=priority[seq1_type],
                                  rank2=priority[seq2_type])
dt.L.clean <- dt.L.clean%>%mutate(pair_type=case_when(
  rank1<rank2~seq1_type, rank2<rank1~seq2_type, TRUE~seq1_type))

table(dt.H.clean$pair_type)
# Dispersed       Proximal       Tandem           WGD / segmental 
# 5954            3874           15378            9714 
table(dt.L.clean$pair_type)
# Dispersed       Proximal       Singleton        Tandem WGD / segmental 
# 34606           30805           58984           91559           51135



dup.colors <- c("WGD / segmental"="#cb7a87", "Tandem"="#e7ba51",
                "Proximal"="#78b43e", "Dispersed"="#999999")

ggplot(dt.H.clean, aes(x=dS, color=pair_type)) + # fill=pair_type, 
  geom_density(linewidth=0.7, alpha=0.9) + 
  # scale_fill_manual(values=dup.colors, name="Duplication mode") +
  scale_color_manual(values=dup.colors, name="Duplication mode") +
  scale_y_continuous(expand=expansion(mult=c(0.00, 0.05))) + 
  scale_x_continuous(expand=expansion(mult=c(0.00, 0.01))) + 
  labs(title=bquote(bolditalic("Oryza sativa")~bold("(H)")), x="dS", y="Density") + 
  theme_custom()