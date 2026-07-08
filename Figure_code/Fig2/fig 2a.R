


library(ggplot2)
library(dplyr)

data <- read.table("merged.txt",
                   header    = FALSE,
                   sep       = " ",
                   col.names = c("CellType", "TotalReadCount"))

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

label_dict <- list(
  "BMem" = "BMem",
  "BNav" = "BNav",
  "CD4EM" = "CD4EM",
  "CD4Naive" = "CD4Nav",
  "CD4Treg" = "CD4Treg",
  "CD8GZMH" = "CD8GZMH",
  "CD8GZMK" = "CD8GZMK",
  "CD8Naive" = "CD8Naive",
  "MAIT" = "MAIT",
  "NKBright" = "NKBright",
  "NKDim" = "NKDim",
  "MonocM" = "MonocM",
  "NoClaM" = "NoClaM"
)

new_celltype_pal <- setNames(
  color_mapping[names(label_dict)],
  unlist(label_dict)
)

data_filtered <- subset(data, CellType %in% names(color_mapping))
data_filtered$CellType <- as.character(sapply(data_filtered$CellType,
                                               function(x) label_dict[[x]]))

median_order <- data_filtered %>%
  group_by(CellType) %>%
  summarise(median_val = median(TotalReadCount), .groups = "drop") %>%
  arrange(desc(median_val)) %>%
  pull(CellType)

data_filtered$CellType <- factor(data_filtered$CellType, levels = median_order)

sample_counts <- data_filtered %>%
  group_by(CellType) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(match(CellType, median_order))

new_labels <- setNames(
  paste0(sample_counts$CellType, "\n(n=", sample_counts$n, ")"),
  sample_counts$CellType
)

overall_median <- median(data_filtered$TotalReadCount)
max_y          <- ceiling(max(data_filtered$TotalReadCount) / 10000) * 10000

p <- ggplot(data_filtered,
            aes(x = CellType, y = TotalReadCount, fill = CellType)) +
  geom_boxplot(
    width         = 0.7,
    outlier.shape = 16,
    outlier.size  = 1.2,
    outlier.alpha = 0.6,
    linewidth     = 0.6
  ) +
  geom_hline(
    yintercept = overall_median,
    color      = "red",
    linetype   = "dashed",
    linewidth  = 0.8
  ) +
  scale_fill_manual(values = new_celltype_pal) +
  scale_x_discrete(labels = new_labels) +
  scale_y_continuous(
    breaks = seq(0, max_y, by = 10000),
    limits = c(0, max_y),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme_minimal() +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y        = element_text(size = 9),
    axis.title.y       = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(color = "black"),
    axis.ticks         = element_line(color = "black"),
    legend.position    = "none"
  ) +
  labs(y = "Total read count", x = NULL)

print(p)

ggsave("total_read_boxplot_sorted.pdf", plot = p, width = 12, height = 5)
ggsave("total_read_boxplot_sorted.png", plot = p, width = 12, height = 5, dpi = 300)

