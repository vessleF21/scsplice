
data <- data.frame(
  cell_type = c("BMem", "BNav", "CD4EM", "CD4Naive", "CD4Treg", "CD8GZMH", 
                "CD8GZMK", "CD8Naive", "MAIT", "NKDim", "NKBright", "MonocM", "NoClaM"),
  overlap = c(58.12, 57.02, 65.57, 62.92, 41.51, 58.23, 54.55, 48.31, 46.43, 60.82, 46.88, 55.08, 54.22),
  sqtl = c(41.88, 42.98, 34.43, 37.08, 58.49, 41.77, 45.45, 51.69, 53.57, 39.18, 53.12, 44.92, 45.78),
  overlap_count = c(68, 69, 400, 431, 22, 138, 84, 57, 13, 222, 15, 65, 45),
  sqtl_count = c(49, 52, 210, 254, 31, 99, 70, 61, 15, 143, 17, 53, 38),
  total = c(117, 121, 610, 685, 53, 237, 154, 118, 28, 365, 32, 118, 83)
)

cell_order <- c("BMem", "BNav", "CD4EM", "CD4Naive", "CD4Treg", "CD8GZMH", 
                "CD8GZMK", "CD8Naive", "MAIT", "NKDim", "NKBright", "MonocM", "NoClaM")

library(ggplot2)
library(dplyr)
library(tidyr)

data_long <- data %>%
  pivot_longer(cols = c(overlap, sqtl), 
               names_to = "type", 
               values_to = "percentage") %>%
  mutate(
    cell_type = factor(cell_type, levels = cell_order),  
    type = factor(type, levels = c("sqtl", "overlap")) 
  )

p1 <- ggplot(data_long, aes(x = cell_type, y = percentage, fill = type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c("sqtl" = "#F28E2B", "overlap" = "#4E79A7"),
                    labels = c("sqtl" = "sQTL", "overlap" = "Overlap"),
                    breaks = c("overlap", "sqtl")) + 
  labs(title = "Proportion of Overlap and sQTL Events by Cell Type",
       x = "Cell Type",
       y = "Percentage (%)",
       fill = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  coord_cartesian(ylim = c(0, 100)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", alpha = 0.5)

print(p1)

data_with_n <- data %>%
  mutate(cell_type_with_n = paste0(cell_type, " (n=", total, ")"))

data_long2 <- data_with_n %>%
  pivot_longer(cols = c(overlap, sqtl), 
               names_to = "type", 
               values_to = "percentage") %>%
  mutate(
    cell_type_with_n = factor(cell_type_with_n, 
                              levels = paste0(cell_order, " (n=", data$total, ")")),
    type = factor(type, levels = c("sqtl", "overlap"))
  )

p2 <- ggplot(data_long2, aes(x = cell_type_with_n, y = percentage, fill = type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = c("sqtl" = "#FC8D62", "overlap" = "#66C2A5"),
                    labels = c("sqtl" = "sQTL", "overlap" = "Overlap"),
                    breaks = c("overlap", "sqtl")) +
  labs(title = "Distribution of Overlap and sQTL Events Across Cell Types",
       x = "Cell Type",
       y = "Percentage (%)",
       fill = "Event Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "top") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50", alpha = 0.5)

