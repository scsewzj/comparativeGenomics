library(tidyr)


array.files <- list.files("mcl/TAGs", pattern = "_arrays.tsv$", full.names = TRUE)

tags.results <- lapply(array.files, function(arr.file) {
  
  # infer spacer range name
  tag.class <- sub("_arrays.tsv$", "", basename(arr.file))
  pair.file <- file.path("mcl/TAGs", paste0(tag.class, "_pairs.tsv"))
  
  # read files
  arrays <- read.delim(arr.file, stringsAsFactors = FALSE)
  pairs  <- read.delim(pair.file, stringsAsFactors = FALSE)
  
  sizes <- as.numeric(arrays$array_size)
  
  data.frame(
    spacer.range = tag.class,
    n.pairs  = nrow(pairs),
    n.arrays = nrow(arrays),
    # size_min  = min(sizes),
    size_mean = mean(sizes),
    size_max  = max(sizes),
    size_sd   = sd(sizes),
    sizes = sizes
  )
})

tags.summary <- do.call(rbind, tags.results)
write.table(tags.summary, file="mcl/TAGs/tags.summary.txt", quote = F, row.names = F, sep="\t")



## ----- Distribution of TAGs sizes

tag.sizes.low <- bind_rows(tags.results[4:6]) %>%
  select(spacer.range, sizes) %>%
  unnest(sizes)
tag.sizes.high <- bind_rows(tags.results[1:3]) %>%
  select(spacer.range, sizes) %>%
  unnest(sizes)
tag.sizes.low$spacer.range <- sub("_TAG_", "/", tag.sizes.low$spacer.range)
tag.sizes.high$spacer.range <- sub("_TAG_", "/", tag.sizes.high$spacer.range)
tag.sizes.low$spacer.range <- sub("_", "-", tag.sizes.low$spacer.range)
tag.sizes.high$spacer.range <- sub("_", "-", tag.sizes.high$spacer.range)

ggplot(tag.sizes.low, aes(x = sizes)) +
  geom_histogram(binwidth = 1, fill = "#6caed6", color = "black") +
  #scale_y_break(c(500, 2000), scales = 0.5, ticklabels = c(500, 2020)) +
  facet_wrap(~ spacer.range, scales = "free_y") +
  labs(x = "TAG array size",
       # title="L (all spacer ranges)",
       y = "# TAG arrays") +
  theme_custom() +
  theme(strip.text = element_text(size = 11, color = "black", face = "bold"))
ggplot(tag.sizes.high, aes(x = sizes)) +
  geom_histogram(binwidth = 1, fill = "#6caed6", color = "black") +
  #scale_y_break(c(300, 1400), scales = 0.5, ticklabels = c(300, 1405)) +
  facet_wrap(~ spacer.range, scales = "free_y") +
  labs(x = "TAG array size",
       #title="H (all spacer ranges)",
       y = "# TAG arrays") +
  theme_custom() +
  theme(strip.text = element_text(size = 11, color = "black", face = "bold"))



tags.size.stats.global <- data.frame(
  Cluster_size = c("total number of arrays", "max size", "2", "3–5", "6-9", ">=10"),
  High = c(length(tag.sizes.high$sizes), max(tag.sizes.high$sizes),
           sum(tag.sizes.high$sizes == 2),
           sum(tag.sizes.high$sizes >= 3 & tag.sizes.high$sizes <= 5),
           sum(tag.sizes.high$sizes >= 6 & tag.sizes.high$sizes <= 9),
           sum(tag.sizes.high$sizes >= 10)),
  Low = c(length(tag.sizes.low$sizes), max(tag.sizes.low$sizes),
          sum(tag.sizes.low$sizes == 2),
          sum(tag.sizes.low$sizes >= 3 & tag.sizes.low$sizes <= 5),
          sum(tag.sizes.low$sizes >= 6 & tag.sizes.low$sizes <= 9),
          sum(tag.sizes.low$sizes >= 10)),
  row.names = NULL)

summarise_tag_sizes <- function(df) {
  df%>%
    group_by(spacer.range)%>%
    summarise(total_arrays = n(),
              max_size = max(sizes),
              size_2 = sum(sizes == 2),
              size_3_5 = sum(sizes >= 3 & sizes <= 5),
              size_6_9 = sum(sizes >= 6 & sizes <= 9),
              size_ge_10 = sum(sizes >= 10),
              .groups = "drop")
}

tags.size.stats.high <- summarise_tag_sizes(tag.sizes.high)
tags.size.stats.low  <- summarise_tag_sizes(tag.sizes.low)


# ---- GO sharing?

list.files("mcl/pantherdb/output/")
go.paths <- c("mcl/pantherdb/output/Slim-low-singletons.txt",
              "mcl/pantherdb/output/Slim-low_TAGs_0.txt",
              "mcl/pantherdb/output/Slim-low_all_TAGs.txt",
              "mcl/pantherdb/output/Slim-low_all_non_TAGs.txt",
              "mcl/pantherdb/output/Slim-high-singletons.txt",
              "mcl/pantherdb/output/Slim-high_TAGs_0.txt",
              "mcl/pantherdb/output/Slim-high_all_TAGs.txt",
              "mcl/pantherdb/output/Slim-high_all_non_TAGs.txt")

tags.gos <- setNames(
  lapply(go.paths, function(go.file) {
    load_go_results(go.file, fdr_filter = 0.05) %>%
      filter(plus_minus == "+") %>%
      pull(GOID)
  }),
  sub("\\.txt$", "", basename(go.paths))
)
names(tags.gos) <- c("L_singletons", "L_TAG/0", "L_TAGs", "L_nonTAGs",
                     "H_singletons", "H_TAG/0", "H_TAGs", "H_nonTAGs")

library(ggvenn)
ggvenn(tags.gos[c(1,2,4)],
       fill_color = c("grey50", "#3674B5", "#4C763B","#e7ba51"),
       stroke_size = 0.5, set_name_size = 4, auto_scale = F, padding = 0.1,fill_alpha=0.5,
       show_percentage=F) + 
  labs(title="(L)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggvenn(tags.gos[c(5,6,8)],
       fill_color = c("grey50", "#3674B5", "#4C763B","#e7ba51"),
       stroke_size = 0.5, set_name_size = 4, auto_scale = F, padding = 0.1, fill_alpha=0.5,
       show_percentage=F) + 
  labs(title="(H)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# common for H
intersect(tags.gos[6][[1]],tags.gos[8][[1]])
# common for L
intersect(tags.gos[2][[1]],tags.gos[4][[1]])


d1 <- load_go_results("mcl/pantherdb/output/Slim-low_all_non_TAGs.txt", fdr_filter = 0.05)
d1 <- d1%>%filter(GOID %in% intersect(tags.gos[2][[1]],tags.gos[4][[1]]))
d2 <- load_go_results("mcl/pantherdb/output/Slim-low_TAGs_0.txt", fdr_filter = 0.05)
d2 <- d2%>%filter(GOID %in% intersect(tags.gos[2][[1]],tags.gos[4][[1]]))

joined.low <- d1%>%select(GOID, GO_term, fold_enrichment_nonTAGs = fold_enrichment, FDR_d1 = FDR) %>%
  inner_join(d2 %>%select(GOID, fold_enrichment_allTAGs = fold_enrichment, FDR_d2 = FDR),
             by = "GOID")
write.table(joined.low, file="mcl/pantherdb/output/common_low_GO_nonTAGs_allTAGs.csv", sep=",", quote=F, row.names = F)

d1 <- load_go_results("mcl/pantherdb/output/Slim-high_all_non_TAGs.txt", fdr_filter = 0.05)
d1 <- d1%>%filter(GOID %in% intersect(tags.gos[6][[1]],tags.gos[8][[1]]))
d2 <- load_go_results("mcl/pantherdb/output/Slim-high_TAGs_0.txt", fdr_filter = 0.05)
d2 <- d2%>%filter(GOID %in% intersect(tags.gos[6][[1]],tags.gos[8][[1]]))

joined.high <- d1%>%select(GOID, GO_term, fold_enrichment_nonTAGs = fold_enrichment, FDR_d1 = FDR) %>%
  inner_join(d2 %>%select(GOID, fold_enrichment_allTAGs = fold_enrichment, FDR_d2 = FDR),
             by = "GOID")
write.table(joined.high, file="mcl/pantherdb/output/common_high_GO_nonTAGs_allTAGs.csv", sep=",", quote=F, row.names = F)





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



