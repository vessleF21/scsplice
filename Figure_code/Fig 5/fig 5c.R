rm(list = ls())
gc()
library(ggplot2)
library(scatterpie)
library(reshape2)
library(dplyr)
library(tidyr)

auto_3 <- read.table("auto_3", header = FALSE, sep = " ", stringsAsFactors = FALSE)
colnames(auto_3) <- c("gene", "celltype", "disease")

all_diseases <- c("CD", "GD", "MS", "RA", "SLE", "T1D", "UC")
x_order <- c("BNav", "BMem", "CD4Naive", "CD4EM", "CD4Treg", "CD8Naive", 
             "CD8GZMH", "CD8GZMK", "MAIT", "NKDim", "NKBright", "MonocM", "NoClaM")
y_order <- sort(unique(auto_3$gene))

auto_3 <- auto_3 %>%
  mutate(
    gene_coord = match(gene, y_order),
    celltype_coord = match(celltype, x_order)
  )

missing_celltypes <- setdiff(unique(auto_3$celltype), x_order)

auto_df_wide <- auto_3 %>%
  group_by(gene_coord, celltype_coord, disease) %>%
  summarise(count = n(), .groups = "drop") %>%
  complete(
    gene_coord = 1:length(y_order),
    celltype_coord = 1:length(x_order),
    disease = all_diseases,
    fill = list(count = 0)
  ) %>%
  pivot_wider(
    names_from = disease,
    values_from = count,
    values_fill = 0
  ) %>%
  select(gene_coord, celltype_coord, all_of(all_diseases)) %>%
  filter(!is.na(gene_coord), !is.na(celltype_coord))

disease_col <- c(
  "CD" = "#fc8d62", "GD" = "#8da0cb", "MS" = "#5CACEE",
  "RA" = "#e78ac3", "SLE" = "#a6d854", "T1D" = "#66c2a5", "UC" = "#FF1493"
)

p <- ggplot() +
  geom_scatterpie(
    aes(x = celltype_coord, y = gene_coord, group = interaction(gene_coord, celltype_coord)),
    data = auto_df_wide,
    cols = all_diseases,
    pie_scale = 2,
    color = NA,
    alpha = 0.8
  ) +
  coord_fixed() +
  scale_x_continuous(
    name = "Cell Type",
    breaks = 1:length(x_order),
    labels = x_order,
    expand = c(0.05, 0.05)
  ) +
  scale_y_continuous(
    name = "Gene",
    breaks = 1:length(y_order),
    labels = y_order,
    expand = c(0.05, 0.05)
  ) +
  scale_fill_manual(values = disease_col, name = "Disease") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Gene-Cell Type Associations by Disease")
print(p)
ggsave("Auto_Gene-Cell Type Associations by Disease.pdf", plot = p, width = 10, height = 8, dpi = 300)
