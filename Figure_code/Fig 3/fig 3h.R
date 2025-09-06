rm(list = ls())
gc()
library(ggplot2)

data <- data.frame(
  Category = c("BMem", "BNav", "BNav", "CD4EM", "CD4Naive", "CD8GZMK", "CD8GZMK", "CD8Naive", "NKBright", "NoClaM"),
  Gender = c("Male", "FeMale", "Male", "Male", "FeMale", "FeMale", "Male", "FeMale", "Male", "Male"),
  Count = c(1, 5, 2, 1, 1, 1, 2, 3, 1, 1)
)

summary_data <- aggregate(Count ~ Category + Gender, data, sum)

summary_data$Category <- factor(summary_data$Category, 
                                levels = c("BNav", "BMem", "CD4Naive", "CD4EM", 
                                           "CD8Naive", "CD8GZMK", "NKBright", "NoClaM"))

sex <- ggplot(summary_data, aes(x = Category, y = Count, color = Gender)) +
  geom_point(size = 5) +
  labs(title = "Category Count by Gender (Dot Plot)",
       x = "Category",
       y = "Count") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid = element_blank(),  
        axis.line = element_line(color = "black", size = 1), 
        axis.ticks = element_line(color = "black", size = 1)) +  
  scale_y_continuous(
    expand = expansion(mult = c(0.1, 0.1)),  
    limits = c(0, max(summary_data$Count) + 1),  
    breaks = seq(0, max(summary_data$Count) + 1, by = 1))
ggsave('sex.pdf', plot = sex, width = 10, height = 10)
