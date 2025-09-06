rm(list = ls())
gc()
library(tidyverse)
library(stringr)
library(data.table)
library(dplyr)
library(ggpmisc)

cell_type_summary_1 <- read.table("summary_1.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
cell_type_summary_2 <- read.table("summary_2.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")
donor_list_2 <- read.table("donor.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

color <- read.table('ref_palette.new.txt', comment.char = '')
colors <- color$V2
names(colors) <- color$V1
plot_2 <- read.table("plot_2.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "")

merged_data <- merge(cell_type_summary_1, plot_2, by = "cell_type")

color_mapping <- c(
  "BNav" = "#01579B", "BMem" = "#A5EDFF", "CD4EM" = "#1B5E20", "CD4Naive" = "#8BC34A",
  "CD4Treg" = "#D4E157", "CD8GZMH" = "#FF5722", "CD8GZMK" = "#FFCDD2", "CD8Naive" = "#F9A825",
  "MAIT" = "#FFEE58", "MonocM" = "#6A1B9A", "NoClaM" = "#E1BEE7", "NKBright" = "#82581F", "NKDim" = "#D7CCC8"
)

merged_data$color <- color_mapping[merged_data$cell_type]

colnames(merged_data) <- c("cell_type","n_donors", "mean_cells","median_cells", "total_cells",
                           "mean_libsize", "median_libsize","total_libsize", "sgene","mean_slope",
                           "median_slope","color","n_cis","n_trans","n_donors.y")    

formula_1 <- y ~ x
correlation <- cor(merged_data$n_donors, merged_data$n_trans, method = "pearson")
correlation_label <- paste0("Pearson's r = ", round(correlation, 2))

p1 <- ggplot(merged_data, aes(x = n_donors, y = n_trans)) +
  geom_point(aes(color = cell_type), size = 5, shape = 21, fill = "white", stroke = 2) +
  scale_color_manual(values = merged_data$color) +
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

ggsave("transgene_vs_donor.pdf", p1, height = 7, width = 7)

formula_2 <- y ~ x
correlation <- cor(merged_data$mean_libsize, merged_data$n_trans, method = "pearson")
correlation_label <- paste0("Pearson's r = ", round(correlation, 2))

p2 <- ggplot(merged_data, aes(x = mean_libsize, y = n_trans)) +
  geom_point(aes(color = cell_type), size = 5, shape = 21, fill = "white", stroke = 2) +
  scale_color_manual(values = merged_data$color) +
  scale_y_log10() +
  stat_poly_eq(formula = y ~ x, aes(label = ..eq.label..), size = 4, label.x = 3, label.y = 0.1, color = "darkred") +
  stat_poly_line(formula = y ~ x, color = "darkblue", size = 1.5) +
  annotate("text", x = Inf, y = Inf, label = correlation_label, hjust = 1.1, vjust = 1.5, size = 5, color = "darkgreen") +
  theme_classic() +
  xlab("Junction read counts") +
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

ggsave("transgene_vs_mean_libsize.pdf", p2, height = 7, width = 7)
