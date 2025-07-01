library(tidyverse)
library(data.table)
library(fs)
library(stringr)
library(dplyr)
library(ggpmisc)
dat_1 <- read.table("dat_1.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")


dat_2 <- dat_1 %>%
  mutate(sqtl = ifelse(sqtl >= 3, ">=3", as.character(sqtl))) %>%
  mutate(sqtl = factor(sqtl, levels = c("1", "2", ">=3")))

dat_3 <- dat_2 %>%
  group_by(cell_type) %>%
  summarize(sgene = n())


cell_type_order <- dat_3 %>%
  arrange(desc(sgene)) %>%
  pull(cell_type)


custom_sqtl_color <- c(
  "1" = "#FDE725", "2" = "#7AD151",
  ">3" = "#22A884"
)    
norm_const <- max(dat_3$sgene) * 1.10
out_dir <- "sum12/"
dir.create(out_dir, showWarnings = FALSE)
setwd("sum12/")

cell_type_order <- c("BNav", "BMem", "CD4Naive", "CD4EM", "CD4Treg", 
                     "CD8Naive", "CD8GZMH", "CD8GZMK", "MAIT", 
                     "NKDim", "NKBright", "MonocM", "NoClaM")

prop_data <- dat_2 %>%
  group_by(cell_type, sqtl) %>%
  summarise(count = n(), .groups = "drop") %>%
  complete(cell_type, sqtl, fill = list(count = 0)) %>%
  left_join(dat_3, by = "cell_type") %>%
  mutate(
    proportion = count / sgene,
    cell_type = factor(cell_type, levels = cell_type_order),
    sqtl = factor(sqtl, levels = c("1", "2", ">=3"))
  )

p <- ggplot() +
  geom_bar(
    data = prop_data,
    aes(x = cell_type, y = proportion, fill = sqtl),
    stat = "identity", 
    width = 0.7, 
    color = "white", 
    size = 0.3
  ) +
  geom_point(
    data = dat_3,
    aes(x = cell_type, y = 1.05, size = sgene),
    color = "#6A3D9A",
    shape = 18,
    show.legend = FALSE
  ) +
  geom_text(
    data = dat_3,
    aes(x = cell_type, y = 1.05, label = sgene),
    color = "#6A3D9A",
    vjust = -1.2,
    size = 3.5,
    fontface = "bold"
  ) +
  coord_polar(theta = "x") +
  ylim(-0.2, 1.3) +
  scale_fill_manual(
    values = c("1" = "#FDE725", "2" = "#7AD151", ">=3" = "#22A884"),
    name = "sQTL/sGene",
    labels = c("1", "2", "≥3")
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
    label = "Purple diamonds show total sGenes per cell type",
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
  ) +
  labs(
    title = "sQTL Distribution Across Cell Types (Circular Stacked Bar Chart)",
    subtitle = "Segments show proportional distribution, diamond size indicates total sGene count",
    caption = "Data Source: sQTL Analysis Pipeline | Outer circle = 100% proportion"
  )


ggsave(
  filename = "sqtl_circular_stacked_bar.pdf",
  plot = p,
  width = 8,  
  height = 6,  
  units = "in",  
  dpi = 300 
)
