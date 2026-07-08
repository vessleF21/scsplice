


rm(list = ls())
gc()

library(ggplot2)
library(dplyr)
library(data.table)

INPUT_FILE <- "auto_3"
OUT_FILE   <- "auto3_bubble.pdf"

CELL_ORDER <- c(
  "BNav", "BMem", "CD4Naive", "CD4EM", "CD4Treg",
  "CD8Naive", "CD8GZMH", "CD8GZMK", "MAIT",
  "NKDim", "NKBright", "MonocM", "NoClaM"
)

DISEASE_ORDER <- c("RA", "CD", "UC", "GD", "T1D", "SLE", "MS")

COLORS_DISEASE <- c(
  "RA"  = "#e41a1c",
  "CD"  = "#377eb8",
  "UC"  = "#4daf4a",
  "GD"  = "#984ea3",
  "T1D" = "#ff7f00",
  "SLE" = "#a65628",
  "MS"  = "#f781bf"
)

df <- fread(INPUT_FILE, header = FALSE,
            col.names = c("gene", "cell_type", "disease"))

bubble_df <- df %>%
  group_by(gene, cell_type) %>%
  summarise(
    n_disease      = n_distinct(disease),
    diseases_label = paste(sort(unique(disease)), collapse = "/"),
    .groups = "drop"
  )

gene_order <- df %>%
  count(gene, sort = TRUE) %>%
  pull(gene)

bubble_df$gene      <- factor(bubble_df$gene,      levels = rev(gene_order))
bubble_df$cell_type <- factor(bubble_df$cell_type, levels = CELL_ORDER)

bubble_df <- bubble_df %>%
  left_join(df %>% distinct(gene, cell_type, disease), by = c("gene", "cell_type")) %>%
  mutate(dot_color = ifelse(n_disease == 1, disease, "Multiple"))

COLORS_DOT <- c(COLORS_DISEASE, "Multiple" = "#888888")

p_bubble <- ggplot(bubble_df, aes(x = cell_type, y = gene)) +
  geom_point(
    aes(size = n_disease, fill = dot_color),
    shape = 21, color = "white", stroke = 0.4
  ) +
  scale_fill_manual(values = COLORS_DOT, name = "Disease") +
  scale_size_continuous(
    range  = c(3, 9),
    breaks = 1:4,
    name   = "No. of diseases"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y      = element_text(face = "bold.italic", color = "black"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  ) +
  labs(x = NULL, y = NULL)

ggsave(OUT_FILE, p_bubble, width = 8, height = 7)