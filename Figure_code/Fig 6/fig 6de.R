ggplot(df, aes(x = genotype, y = level)) +
  geom_violin(
    aes(group = genotype),
    fill = NA,
    color = "darkgreen",
    draw_quantiles = c(0.5)
  ) +
  geom_boxplot(
    width = 0.1,
    fill = "white",
    color = "darkgreen",
    outlier.shape = NA
  ) +
  stat_summary(
    fun = mean,
    geom = "line",
    aes(group = 1),
    color = "red",
    size = 1
  ) +
  annotate(
    "text",
    x = 1,
    y = 4,
    label = paste0("P = ", format(p_val, scientific = TRUE, digits = 2)),
    hjust = 0
  ) +
  labs(
    x = "snp",
    y = "intron"
  ) +
  theme_classic()
