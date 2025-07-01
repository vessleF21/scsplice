rm(list = ls())
gc()
library(tidyverse)
library(stringr)
library(data.table)
library(dplyr)
library(ggpmisc)

cell_type_summary_1 <- read.table("summary_1.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
cell_type_summary_2 <- read.table("summary_2.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
donor_list_2 <- read.table("donor_list.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

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

cell_type_summary_1$color <- color_mapping[cell_type_summary_1$cell_type]

formula_1 <- y ~ x
correlation <- cor(cell_type_summary_1$n_donors, cell_type_summary_1$sgene, method = "pearson")
correlation_label <- paste0("Pearson's r = ", round(correlation, 2))

p1 <- ggplot(cell_type_summary_1, aes(x = n_donors, y = sgene)) +
  geom_point(aes(color = cell_type), size = 5, shape = 21, fill = "white", stroke = 2) +
  scale_color_manual(values = cell_type_summary_1$color) +
  scale_y_log10() +
  stat_poly_eq(formula = y ~ x, aes(label = ..eq.label..), size = 4, label.x = 3, label.y = 0.1, color = "darkred") +
  stat_poly_line(formula = y ~ x, color = "darkblue", size = 1.5) +
  annotate("text", x = Inf, y = Inf, label = correlation_label, hjust = 1.1, vjust = 1.5, size = 5, color = "darkgreen") +
  theme_classic() +
  xlab("Number of donors") +
  ylab("Number of sGenes") +
  guides(color = FALSE) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(face = "bold", size = 16, color = "black"),
    axis.text = element_text(size = 22, color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18, color = "darkblue"),
    plot.subtitle = element_text(size = 22, hjust = 0.5, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.line = element_line(color = "black", size = 0.6)
  )

print(p1)
ggsave("sgene_vs_donor.pdf", p1, height = 7, width = 7)


formula_2 <- y ~ x
correlation <- cor(cell_type_summary_1$mean_cells, cell_type_summary_1$sgene, method = "pearson")
correlation_label <- paste0("Pearson's r = ", round(correlation, 2))

p2 <- ggplot(cell_type_summary_1, aes(x = mean_cells, y = sgene)) +
  geom_point(aes(color = cell_type), size = 5, shape = 21, fill = "white", stroke = 2) +
  scale_color_manual(values = cell_type_summary_1$color) +
  stat_poly_eq(formula = formula_2, aes(label = ..eq.label..), size = 4, label.x = 3, label.y = 0.1, color = "darkred") +
  stat_poly_line(formula = formula_2, color = "darkblue", size = 1.5) +
  annotate("text", x = Inf, y = Inf, label = correlation_label, hjust = 1.1, vjust = 1.5, size = 5, color = "darkgreen") +
  theme_classic() +
  xlab("Cells per donor") +
  ylab("Number of sGenes") +
  guides(color = FALSE) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(face = "bold", size = 16, color = "black"),
    axis.text = element_text(size = 22, color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18, color = "darkblue"),
    plot.subtitle = element_text(size = 22, hjust = 0.5, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.line = element_line(color = "black", size = 0.6)
  )

print(p2)
ggsave("sgene_vs_mean_cells.pdf", p2, height = 7, width = 7)


correlation <- cor(cell_type_summary_1$n_donors, cell_type_summary_1$mean_slope, method = "pearson")
correlation_label <- paste0("Pearson's r = ", round(correlation, 2))
formula_5 <- y ~ x + poly(x, 2)

p3 <- ggplot(cell_type_summary_1, aes(x = n_donors, y = mean_slope)) +
  geom_point(aes(color = cell_type), size = 5, shape = 21, fill = "white", stroke = 2) +
  scale_color_manual(values = cell_type_summary_1$color) +
  xlab("Number of donors") +
  ylab("Mean absolute effect size") +
  stat_poly_eq(formula = formula_5, aes(label = ..eq.label..), label.x = "right", size = 4, color = "darkred") +
  stat_poly_line(formula = formula_5, color = "darkblue", size = 1.5) +
  annotate("text", x = Inf, y = Inf, label = correlation_label, hjust = 1.1, vjust = 1.5, size = 5, color = "darkgreen") +
  theme_classic() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(face = "bold", size = 16, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18, color = "darkblue"),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.line = element_line(color = "black", size = 0.6)
  )

print(p3)
ggsave("effect_size_vssample_size_with_labels.pdf", p3, height = 7, width = 7)

correlation <- cor(cell_type_summary_1$mean_libsize, cell_type_summary_1$mean_slope, method = "pearson")
correlation_label <- paste0("Pearson's r = ", round(correlation, 2))
formula_6 <- y ~ x

p4 <- ggplot(cell_type_summary_1, aes(x = mean_libsize, y = mean_slope)) +
  geom_point(aes(color = cell_type), size = 5, shape = 21, fill = "white", stroke = 2) +
  scale_color_manual(values = cell_type_summary_1$color) +
  xlab("Library size") +
  ylab("Mean absolute effect size") +
  stat_poly_eq(formula = formula_6, aes(label = ..eq.label..), label.x = "right", size = 4, color = "darkred") +
  stat_poly_line(formula = formula_6, color = "darkblue", size = 1.5) +
  annotate("text", x = Inf, y = Inf, label = paste0("Pearson's r = ", round(correlation, 2)), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "darkgreen") +
  theme_classic() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(face = "bold", size = 16, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18, color = "darkblue"),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.line = element_line(color = "black", size = 0.6)
  )

print(p4)
ggsave("effect_size_vs_library_size.pdf", p4, height = 7, width = 7)
