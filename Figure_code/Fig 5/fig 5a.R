library(ComplexHeatmap)
celltypes <- read.table('cell_type.txt')$V1
color <- read.table('palette.txt', comment.char = '')
palatte <- color$V2
names(palatte) <- color$V1


celltype_pal <- c(
  "BNav" = "#01579B",
  "BMem" = "#A5EDFF",
  "CD4EM" = "#1B5E20",
  "CD4Naive" = "#8BC34A",
  "CD4Treg" = "#D4E157",
  "CD8GZMH" = "#FF5722",
  "CD8GZMK" = "#FFCDD2",
  "CD8Naive" = "#F9A825",
  "MAIT" = "#FFEE58",
  "MonocM" = "#6A1B9A",
  "NoClaM" = "#E1BEE7",
  "NKBright" = "#82581F",
  "NKDim" = "#D7CCC8"
)  

label_dict<-list("BMem"="BMem",
                 "BNav"="BNav",
                 "CD4EM"="CD4EM",
                 "CD4Naive"="CD4Naive",
                 "CD4Treg"="CD4Treg",
                 "CD8GZMH"="CD8GZMH",
                 "CD8GZMK"="CD8GZMK",
                 "CD8Naive"="CD8Naive",
                 "MAIT"="MAIT",
                 "MonocM"="MonocM",
                 "NKBright"="NKBright",
                 "NKDim"="NKDim",
                 "NoClaM"="NoClaM")

dat <- read.delim('trans_sQTL.tsv')


dat_subset <- dat[, c("gene_name", "cell_type")]


library(dplyr)

gene_cell_summary <- dat_subset %>%
  group_by(gene_name) %>%
  summarise(n_cell_types = n_distinct(cell_type))

dat_tagged <- dat_subset %>%
  left_join(gene_cell_summary, by = "gene_name") %>%
  mutate(dup_flag = ifelse(n_cell_types > 1, "rep", "norep"))

dup_stats <- dat_tagged %>%
  group_by(cell_type, dup_flag) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cell_type) %>%
  mutate(
    total = sum(count),
    proportion = count / total
  )

dup_stats
dup_stats <- dup_stats %>%
  mutate(dup_flag = factor(dup_flag, levels = c("norep", "rep")))
dup_stats <- dup_stats %>%
  mutate(dup_flag = factor(dup_flag, levels = c("norep", "rep")))

p <- ggplot() +
  geom_bar(
    data = dup_stats,
    aes(x = cell_type, y = proportion, fill = dup_flag),
    stat = "identity", 
    position = "stack", 
    width = 0.7, 
    color = "white", 
    size = 0.3
  ) +

  geom_point(
    data = dup_stats %>% distinct(cell_type, .keep_all = TRUE),
    aes(x = cell_type, y = 1.05, size = total),
    color = "#6A3D9A",
    shape = 18,
    show.legend = FALSE
  ) +

  geom_text(
    data = dup_stats %>% distinct(cell_type, .keep_all = TRUE),
    aes(x = cell_type, y = 1.05, label = total),
    color = "#6A3D9A",
    vjust = -1.2,
    size = 3.5,
    fontface = "bold"
  ) +

  coord_polar(theta = "x") +

  ylim(-0.2, 1.3) +

  scale_fill_manual(
    values = c("norep" = "#A6D854", "rep" = "#1B9E77"), 
    name = "sQTL",
    labels = c("norep (1 sQTL)", "rep (≥2 sQTLs)"),
    drop = FALSE
  ) +

  scale_size_continuous(range = c(3, 8)) +

  geom_hline(
    yintercept = c(0.25, 0.5, 0.75, 1), 
    color = "gray90", 
    linetype = "dashed", 
    size = 0.3
  ) +

  annotate(
    "text", 
    x = 0.5, 
    y = c(0.25, 0.5, 0.75, 1), 
    label = c("25%", "50%", "75%", "100%"),
    color = "gray40", 
    size = 3.5
  ) +

  annotate(
    "text", 
    x = 0.5, 
    y = 0.1, 
    label = "sGenes",
    color = "#6A3D9A", 
    size = 3.5
  ) +

  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, color = "gray30", margin = margin(b = 10)),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 10, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    plot.caption = element_text(color = "gray40", hjust = 0.5, margin = margin(t = 10)),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )
  )

print(p)
ggsave("sQTL_distribution_circular.png", width = 10, height = 10, dpi = 300)

