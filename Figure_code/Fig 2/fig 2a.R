
library(ggplot2)
library(dplyr)

data <- read.table("merged.txt", header = FALSE, sep = " ", col.names = c("CellType", "TotleReadCount"))

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

data_filtered <- data_filtered %>%
  group_by(CellType) %>%
  mutate(MedianTotleReadCount = median(TotleReadCount)) %>%
  ungroup()

data_filtered$CellType <- factor(data_filtered$CellType, 
                                 levels = data_filtered %>%
                                   arrange(desc(MedianTotleReadCount)) %>%
                                   pull(CellType) %>%
                                   unique())

levels(data_filtered$CellType) <- sapply(levels(data_filtered$CellType), function(x) label_dict[[x]])

p <- ggplot(data_filtered, aes(x = CellType, y = TotleReadCount, fill = CellType)) +
  geom_boxplot(color = "black") +
  scale_fill_manual(values = new_celltype_pal) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  labs(
    title = "Totle read count by Cell Type", 
    y = "Totle read count", 
    x = "Cell Type"
  )

ggsave("Totle_read_count.pdf", plot = p, width = 10, height = 6, dpi = 300)






