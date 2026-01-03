

# ----- load the pairs
tags.pairs <- setNames(lapply(array.files, function(arr.file) {
  tag.class <- sub("_arrays.tsv$", "", basename(arr.file))
  pair.file <- file.path("mcl/TAGs", paste0(tag.class, "_pairs.tsv"))
  read.delim(pair.file, stringsAsFactors = FALSE)}),
  sub("_arrays.tsv$", "", basename(array.files))
)
names(tags.pairs)
head(tags.pairs[[1]])
head(dt.H.clean[,1:3])
head(dt.L.clean[,1:3])



ds.H <- dt.H.clean %>%
  transmute(g1 = pmin(seq1, seq2),
            g2 = pmax(seq1, seq2),
            dS = dS)

ds.L <- dt.L.clean %>%
  transmute(g1 = pmin(seq1, seq2),
            g2 = pmax(seq1, seq2),
            dS = dS)


tags.pairs <- lapply(names(tags.pairs), function(name) {
  ds_lookup <- if (startsWith(name, "H_")) ds.H else ds.L
  tags.pairs[[name]] %>%
    mutate(
      g1 = pmin(gene1, gene2),
      g2 = pmax(gene1, gene2)
    ) %>%
    left_join(ds_lookup, by = c("g1", "g2")) %>%
    select(-g1, -g2)
}) %>% setNames(names(tags.pairs))

# Check for missing dS
sapply(names(tags.pairs), function(name) {
  n_missing <- sum(is.na(tags.pairs[[name]]$dS))
  cat(sprintf("%s: %d/%d pairs missing dS\n", 
              name, n_missing, nrow(tags.pairs[[name]])))
  n_missing
})



ggplot() + 
  geom_density(data=dt.H.clean, mapping=aes(x=dS, color="all duplicated pairs"), 
               linewidth=0.7, alpha=0.7) + 
  geom_density(data=tags.pairs[["H_TAG_0"]], mapping=aes(x=dS, color="TAG:0 pairs"), 
               linewidth=0.7, alpha=0.7) + 
  geom_density(data=tags.pairs[["H_TAG_1_5"]], mapping=aes(x=dS, color="TAG:1-5 pairs"), 
               linewidth=0.7, alpha=0.7) + 
  geom_density(data=tags.pairs[["H_TAG_6_10"]], mapping=aes(x=dS, color="TAG:6-10 pairs"), 
               linewidth=0.7, alpha=0.7) + 
  geom_density(data=nonTAGs.H, mapping=aes(x=dS, color="non-TAG pairs"), 
               linewidth=0.7, alpha=0.7) + 
  scale_color_manual(name=" ", values=c("all duplicated pairs"="black", "TAG:0 pairs"="#4C763B",
                                        "TAG:1-5 pairs"="#1C6EA4", "TAG:6-10 pairs"="#e7ba51",
                                        "non-TAG pairs"="grey")) +
  scale_y_continuous(expand=expansion(mult=c(0.00, 0.05))) + 
  scale_x_continuous(expand=expansion(mult=c(0.00, 0.01))) + 
  labs(title=bquote(bolditalic("Oryza sativa")~bold("(H)")), x="dS", y="Density") + 
  theme_custom()

ggplot() + 
  geom_density(data=dt.L.clean, mapping=aes(x=dS, color="all duplicated pairs"), 
               linewidth=0.7, alpha=0.7) + 
  geom_density(data=tags.pairs[["L_TAG_0"]], mapping=aes(x=dS, color="TAG:0 pairs"), 
               linewidth=0.7, alpha=0.7) + 
  geom_density(data=tags.pairs[["L_TAG_1_5"]], mapping=aes(x=dS, color="TAG:1-5 pairs"), 
               linewidth=0.7, alpha=0.7) + 
  geom_density(data=tags.pairs[["L_TAG_6_10"]], mapping=aes(x=dS, color="TAG:6-10 pairs"), 
               linewidth=0.7, alpha=0.7) + 
  geom_density(data=nonTAGs.L, mapping=aes(x=dS, color="non-TAG pairs"), 
               linewidth=0.7, alpha=0.7) + 
  scale_color_manual(name=" ", values=c("all duplicated pairs"="black", "TAG:0 pairs"="#4C763B",
                                        "TAG:1-5 pairs"="#1C6EA4", "TAG:6-10 pairs"="#e7ba51",
                                        "non-TAG pairs"="grey")) +
  scale_y_continuous(expand=expansion(mult=c(0.00, 0.05))) + 
  scale_x_continuous(expand=expansion(mult=c(0.00, 0.01))) + 
  labs(title=bquote(bolditalic("Oryza sativa")~bold("(L)")), x="dS", y="Density") + 
  theme_custom()


# we need non-TAG!
# 1142 non-TAG pairs for H
# 4275 non-TAG pairs for L

nonTAGs.H <- read.csv("mcl/TAGs/H_non_TAG_pairs.tsv", sep="\t")
nonTAGs.L <- read.csv("mcl/TAGs/L_non_TAG_pairs.tsv", sep="\t")


nonTAGs.H <- nonTAGs.H %>%
  mutate(key1 = pmin(gene1, gene2),
         key2 = pmax(gene1, gene2)) %>%
  left_join(ds.H,
            by = c("key1" = "g1", "key2" = "g2")) %>%
  select(-key1, -key2)
sum(is.na(nonTAGs.H$dS)) # 235 out of 1142

nonTAGs.L <- nonTAGs.L %>%
  mutate(key1 = pmin(gene1, gene2),
         key2 = pmax(gene1, gene2)) %>%
  left_join(ds.L,
            by = c("key1" = "g1", "key2" = "g2")) %>%
  select(-key1, -key2)
sum(is.na(nonTAGs.L$dS)) # 1154 out of 4275


