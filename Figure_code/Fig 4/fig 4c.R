library(ggplot2)
library(dplyr)

data <- data.frame(
  tissues = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13),
  count = c(330, 176, 72, 35, 24, 12, 11, 10, 4, 6, 4, 2, 6)
)

data <- data %>% 
  arrange(tissues) %>%
  mutate(
    label = paste0(tissues, " (", count, ")"),
    tissues_factor = factor(tissues, levels = sort(unique(tissues)))
  )

colors <- c(
  "#B3E5E1",
  "#FFD4A3", 
  "#FFDFD4",
  "#E8D4F2",
  "#FFE8C6",
  "#D4E8D4",
  "#C6D4E8",
  "#F2D4D4",
  "#E8E8D4",
  "#D4F2F2",
  "#F2E8D4",
  "#E8D4E8",
  "#D4D4E8"
)

p <- ggplot(data, aes(x = "", y = count, fill = tissues_factor)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 0.5) +
  coord_polar(theta = "y", start = 0) +
  scale_fill_manual(
    values = colors,
    labels = data$label,
    name = "Number of shared tissues"
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  annotate("text", x = 0, y = 0, 
           label = "381,202\noverlapped\ngene pairs",
           size = 4, fontface = "bold")

print(p)
