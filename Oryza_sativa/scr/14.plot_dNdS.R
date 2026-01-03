
setwd("/home/ngoc/comparativeGenomics/Oryza_sativa/")

library(ggplot2)
library(dplyr)


# ----- Loading dS results for the high-stringent dataset
dt.H <- read.csv("paml-workflow/high/dNdS.csv")
dt.H.clean <- dt.H%>%filter(dS<=6)  


# ----- Loading dS results for the low-stringent dataset
dt.L <- read.csv("paml-workflow/low/dNdS.csv")
dt.L.clean <- dt.L%>%filter(dS<6)



## ----- dS and dS_SE?
ggplot(dt.H.clean, aes(x = dS, y = dS_SE)) +
  geom_point(alpha = 0.2) +
  labs(x="dS", y="dS_SE", title="(H) dS<6, 18893 pairs") +
  theme_custom()
ggplot(dt.H.clean%>%filter(dS_SE<5), aes(x = dS, y = dS_SE)) +
  geom_point(alpha = 0.2) +
  labs(x="dS", y="dS_SE", title="(H) dS<6 & SE<5, 13876 pairs") +
  theme_custom()
ggplot(dt.L.clean, aes(x = dS, y = dS_SE)) +
  geom_point(alpha = 0.2) +
  labs(x="dS", y="dS_SE", title="(L) dS<6, 155701 pairs") +
  theme_custom()
ggplot(dt.L.clean%>%filter(dS_SE<5), aes(x = dS, y = dS_SE)) +
  geom_point(alpha = 0.2) +
  labs(x="dS", y="dS_SE", title="(L) dS<6 & SE<5, 73962 pairs") +
  theme_custom()


dt.H.clean <- dt.H%>%filter(dS<=6 & dS_SE<5) 
dt.L.clean <- dt.L%>%filter(dS<6  & dS_SE<5)
rm(dt.L)
rm(dt.H)


## ----- Histogram
ggplot(dt.H.clean%>%filter(dS_SE<10), aes(x=dS)) + 
  geom_histogram(bins=100, fill="#4C763B", alpha=0.7) + # , color="black"
  scale_y_continuous(expand=expansion(mult=c(0.0, 0.1))) + 
  scale_x_continuous(expand=expansion(mult=c(0.0, 0.01))) + 
  # geom_vline(xintercept = c(0.18, 1.3, 1.8, 3.5), linetype = "dotted", size = 0.5,
  #            color = c("red3", "black", "black", "red3")) +
  # geom_text(data = data.frame(x = c(0.18, 1.3, 1.8, 3.5)), aes(x = x, y = 450, label = x),
  #           inherit.aes = FALSE, vjust = -0.3, hjust = -0.3, size = 3.5, color="red3") +
  labs(title=bquote(bolditalic("Oryza sativa")~bold("(H)")), x="dS", y="# duplicated genes") + 
  theme_custom()
ggplot(dt.L.clean, aes(x=dS)) + 
  geom_histogram(bins=120, fill="#6caed6", alpha=0.7) + 
  scale_y_continuous(expand=expansion(mult=c(0, 0.1))) + 
  scale_x_continuous(expand=expansion(mult=c(0, 0.01))) +
  # geom_vline(xintercept = c(1.4, 2.4, 3.8), linetype = "dotted", size = 0.5,
  #            color = c("black", "black", "red3")) +
  # geom_text(data = data.frame(x = c(1.4, 2.4, 3.8)), aes(x = x, y = 5500, label = x),
  #           inherit.aes = FALSE, vjust = -0.3, hjust = -0.3, size = 3.5, color="red3") +
  labs(title=bquote(bolditalic("Oryza sativa")~bold("(L)")), x="dS", y="# duplicated genes") + 
  theme_custom()



## ----- Density curve
ggplot() + 
  geom_density(data=dt.H.clean, mapping=aes(x=dS, fill="High"), 
               color="#4C763B", linewidth=0.5, alpha=0.7) + 
  geom_density(data=dt.L.clean, mapping=aes(x=dS, fill="Low"), 
               color="#6caed6", linewidth=0.5, alpha=0.7) + 
  scale_fill_manual(name=" ", values=c("High"="#4C763B", "Low"="#6caed6")) +
  scale_y_continuous(expand=expansion(mult=c(0.00, 0.05))) + 
  scale_x_continuous(expand=expansion(mult=c(0.00, 0.01))) + 
  labs(title=bquote(bolditalic("Oryza sativa")), x="dS", y="Density") + 
  theme_custom()









## ----- Helper functions
theme_custom <- function() {
  theme_minimal() +
    theme(panel.border= element_blank(), # element_rect(colour="black", fill=NA, linewidth=0.5), # 
          axis.line.x= element_line(colour="black", linewidth=0.5), axis.line.y= element_line(colour="black", linewidth=0.5),
          plot.title=element_text(color="black", face="bold", hjust=0.5, size=13),
          panel.grid.major=element_blank(), #element_line(color="grey97"),
          panel.grid.minor=element_blank(),
          axis.text=element_text(size=11, color="black"),
          axis.text.x=element_text(size=11, color="black", angle=0, hjust=1),
          axis.title=element_text(size=13, color="black", face="bold"),
          legend.title=element_text(size=11, color="black", face="bold"),
          legend.text=element_text(size=11, color="black")) 
}





## ----- Color palettes

c("#B0CE88", "#4C763B", "#043915", "#84994F", "#626F47", "#A4B465", "#D2D0A0",
  "#154D71", "#1C6EA4", "#33A1E0", "#3674B5", "#578FCA",
  "#A66E38", "#F0BB78", "#bd9e39", "#e7ba51", "#e7cb94")

c("#910e24", "#b13d50", "#cb7a87", "#e6b8bf", "#843c3a", "#ac494a", "#d5616c", "#e8969d",
  "#99600f", "#b2823f", "#cbaa7a", "#bd9e39", "#e7ba51", "#e7cb94",
  "#539a10", "#78b43e", "#a3cb7a", "#cfe6b8", "#627a39", "#8ca252", "#b5cf6a", "#cfdb9b",
  "#10829a", "#3f9fb4", "#7abecc", "#b7dee5", "#3180bb", "#6caed6", "#9ecae2", "#c6dcef",
  "#383a77", "#5254a3", "#6c6ece", "#9d9edf", "#756cb1", "#9f9ac7", "#bcbddc", "#dadaea",
  "#e6550b", "#feaf6b", "#fdd0a2", "#30a255", "#74c476", "#a0da9a", "#c6e9bf",
  "#636363", "#969696", "#bdbdbd", "#d9d9d9", "#333333", "#666666", "#999999", "#cccccc")
  