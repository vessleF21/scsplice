rm(list = ls())
gc()

library(tidyverse)
library(data.table)
library(fs)
library(stringr)
library(dplyr)
library(ggpmisc)

multi_sqtl_gene <- read.table("heterogeneity.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
correlation <- cor(multi_sqtl_gene$total_sgene, multi_sqtl_gene$prop_multi_sqtl_gene, method = "pearson")
correlation_label <- paste0("Pearson's r = ", round(correlation, 2))

cell_type_colors <- setNames(multi_sqtl_gene$color, multi_sqtl_gene$cell_type)

color_mapping <- c(
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

p2 <- ggplot(multi_sqtl_gene, aes(x = total_sgene, y = prop_multi_sqtl_gene)) +
  geom_point(aes(color = cell_type), size = 4, shape = 21, fill = "white", stroke = 2) +
  scale_color_manual(values = color_mapping) +
  stat_poly_eq(formula = y ~ x, aes(label = ..eq.label..), size = 4, label.x = 3, label.y = 0.1, color = "darkred") +
  stat_poly_line(formula = y ~ x, color = "darkblue", size = 1.5) +
  annotate("text", x = Inf, y = Inf, label = correlation_label, hjust = 1.1, vjust = 1.5, size = 5, color = "darkgreen") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(face = "bold", size = 30, color = "black"),
    axis.text = element_text(size = 30, color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18, color = "darkblue"),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.line = element_line(color = "black", size = 0.6)
  ) +
  xlab("Number of sGenes") +
  ylab("Proportion of sGenes\nwith >1 sSNPs")

print(p2)
ggsave("heterogeneity.pdf", p2, width = 7, height = 7)
