library(ggplot2)
library(ggpointdensity)

data <- read.table("all_R2.txt", header = FALSE, col.names = c("puQTL", "eQTL"))

cor_test <- cor.test(data$puQTL, data$eQTL, method = "spearman")
cor_value <- round(cor_test$estimate, 2)
p_value <- format.pval(cor_test$p.value, digits = 2, eps = 2.2e-16)

ggplot(data, aes(x = puQTL, y = eQTL)) +
  geom_pointdensity(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE, size = 1) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "puQTL effect size", y = "eQTL effect size") +
  coord_cartesian(xlim = c(-1.6, 1.6), ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  ) +
  annotate("text", x = -1.4, y = 1.8,
           label = bquote(r[s] == .(cor_value)),
           parse = TRUE, size = 5, hjust = 0) +
  annotate("text", x = -1.4, y = 1.6,
           label = paste0("P-value = ", p_value),
           size = 5, hjust = 0)
