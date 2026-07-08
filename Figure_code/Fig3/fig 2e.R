

rm(list = ls())
gc()

library(tidyverse)
library(data.table)
library(ggplot2)
library(scales)

label_dict <- list(
  "BNav"     = "BNav",
  "BMem"     = "BMem",
  "CD4Naive" = "CD4Nav",
  "CD4EM"    = "CD4EM",
  "CD4Treg"  = "CD4Treg",
  "CD8Naive" = "CD8Naive",
  "CD8GZMH"  = "CD8GZMH",
  "CD8GZMK"  = "CD8GZMK",
  "MAIT"     = "MAIT",
  "NKDim"    = "NKDim",
  "NKBright" = "NKBright",
  "MonocM"   = "cM",
  "NoClaM"   = "ncM"
)

color_mapping <- c(
  "BNav"     = "#01579B",
  "BMem"     = "#A5EDFF",
  "CD4EM"    = "#1B5E20",
  "CD4Nav"   = "#8BC34A",
  "CD4Treg"  = "#D4E157",
  "CD8GZMH"  = "#FF5722",
  "CD8GZMK"  = "#FFCDD2",
  "CD8Naive" = "#F9A825",
  "MAIT"     = "#FFEE58",
  "cM"       = "#6A1B9A",
  "ncM"      = "#E1BEE7",
  "NKBright" = "#82581F",
  "NKDim"    = "#D7CCC8"
)

data <- read.csv("q1234", sep = '\t')
dat.filter <- data

for (i in seq_along(label_dict)) {
  dat.filter$Celltype <- gsub(
    paste0("^", names(label_dict)[i], "$"),
    label_dict[[i]],
    dat.filter$Celltype
  )
}

median_order <- dat.filter %>%
  group_by(Celltype) %>%
  summarise(med = median(num), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(Celltype)

dat.filter$Celltype <- factor(dat.filter$Celltype, levels = median_order)

p <- ggplot(dat.filter, aes(x = Celltype, y = num, fill = Celltype)) +
  geom_boxplot(
    width         = 0.7,
    outlier.shape = 16,
    outlier.size  = 1.0,
    outlier.alpha = 0.5,
    linewidth     = 0.5
  ) +
  scale_fill_manual(values = color_mapping) +
  scale_y_continuous(labels = label_comma(), breaks = pretty_breaks(n = 5)) +
  ggtitle("Pseudobulk") +
  theme_classic() +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 13),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y     = element_text(size = 9),
    axis.title.y    = element_text(size = 11),
    legend.position = "none"
  ) +
  labs(y = "Number of AS genes", x = NULL)

print(p)

ggsave("AS_gene_pseudobulk.pdf", plot = p, width = 10, height = 6)
ggsave("AS_gene_pseudobulk.png", plot = p, width = 10, height = 6, dpi = 300)


