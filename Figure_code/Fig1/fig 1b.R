

rm(list = ls())
gc()
library(ggplot2)

cell_data <- read.csv("data1.csv", header = TRUE)

median_val <- median(cell_data$x)
mean_val   <- mean(cell_data$x)
sd_val     <- sd(cell_data$x)

p <- ggplot(cell_data, aes(x = "", y = x)) +
  geom_violin(
    fill      = NA,
    color     = "black",
    linewidth = 0.6
  ) +
  geom_boxplot(
    width         = 0.1,
    fill          = NA,
    color         = "black",
    linewidth     = 0.6,
    outlier.shape = NA
  ) +
  ggplot2::annotate(
    "text",
    x     = 1,
    y     = median_val,
    label = paste0("Median\n", median_val),
    color = "darkred",
    size  = 3.5
  ) +
  labs(
    x = "Sample (n=980)",
    y = "Cell Count"
  ) +
  scale_y_continuous(
    limits = c(0, 3000),
    breaks = c(0, 1000, 2000, 3000),
    expand = c(0, 0)
  ) +
  theme_classic() +
  theme(
    axis.line    = element_line(colour = "black", linewidth = 0.6),
    axis.ticks   = element_line(colour = "black", linewidth = 0.6),
    axis.text.x  = element_text(size = 11, color = "black"),
    axis.text.y  = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12, color = "black",
                                 margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, color = "black",
                                 margin = margin(r = 10)),
    panel.grid   = element_blank()
  )

print(p)

ggsave("cell_count_plot.pdf", p, width = 4, height = 5, dpi = 300)
ggsave("cell_count_plot.png", p, width = 4, height = 5, dpi = 300)
