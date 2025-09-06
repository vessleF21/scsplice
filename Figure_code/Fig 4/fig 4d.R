original_values <- c(1369, 204, 1, 1, 3, 10, 40, 1)
original_percentages <- c(0.00, 100.00, 16.67, 20.00, 25.00, 33.33, 50.00, 66.67)

values <- c(1369, 204, sum(original_values[3:8]))
percentages <- c(0.00, 100.00, sum(original_percentages[3:8]))

labels <- c(
  paste("1369 (", 0.00, "%)", sep=""),
  paste("204 (", 100.00, "%)", sep=""),
  paste("Others (", round(sum(original_percentages[3:8]), 2), "%)", sep="")
)

colors <- c("#FF6B6B", "#4ECDC4", "#45B7D1")

pie(values, 
    labels = labels,
    main = "Pie Chart",
    col = colors,
    cex = 1.2,
    border = "white",
    lwd = 2)

library(ggplot2)

data <- data.frame(
  values = values,
  percentages = percentages,
  labels = labels,
  category = c("Category 1", "Category 2", "Others")
)

ggplot(data, aes(x = "", y = values, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  labs(title = "Optimized Pie Chart") +
  scale_fill_manual(values = c("#FF6B6B", "#4ECDC4", "#45B7D1")) +
  geom_text(aes(label = labels), 
            position = position_stack(vjust = 0.5),
            color = "white",
            fontface = "bold",
            size = 4) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

library(plotly)

plot_ly(data, 
        labels = ~category, 
        values = ~values, 
        type = 'pie',
        textinfo = 'label+percent+value',
        textposition = 'inside',
        marker = list(colors = c("#FF6B6B", "#4ECDC4", "#45B7D1"),
                      line = list(color = '#FFFFFF', width = 2))) %>%
  layout(title = list(text = "Interactive Pie Chart", font = list(size = 18)),
         showlegend = TRUE,
         legend = list(orientation = "v", x = 1.02, y = 0.5))
