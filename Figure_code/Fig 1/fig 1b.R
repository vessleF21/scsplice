rm(list = ls())
gc()
library(ggplot2)
cell_data <- read.csv("data1.csv", header = TRUE)

median_val <- median(cell_data$x)
mean_val <- mean(cell_data$x)
sd_val <- sd(cell_data$x)

p <- ggplot(cell_data, aes(x = "", y = x)) +
  geom_violin(fill = NA, color = "black") +
  geom_boxplot(width = 0.1, 
               fill = NA, 
               color = "black", 
               outlier.shape = NA) +
  annotate("text", 
           x = 1, 
           y = median_val, 
           label = paste0("Median\n", median_val),
           color = "darkred", size = 3.5) +
  labs(title = "Cell Count Distribution",
       x = "",
       y = "Cell Count") +
  theme_minimal()

ggsave("cell_count_plot.pdf", plot = p, width = 5, height = 5)

