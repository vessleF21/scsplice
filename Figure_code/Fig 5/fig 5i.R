ggplot(df, aes(x = cis_logP, y = trans_logP)) +
  geom_point(aes(color = color), size = 2) +
  scale_color_identity() +
  labs(
    x = expression("IVNS1ABP cis-eQTL " * -log[10](italic(P))),
    y = expression("RPS24 trans-sQTL " * -log[10](italic(P)))
  ) +
  theme_classic() +

  geom_text(
    data = subset(df, highlight),
    aes(label = paste0(snp, "\nH4 = 0.99")),
    hjust = -0.1,
    vjust = 1.2,
    size = 4
  ) +
  xlim(0, max(df$cis_logP) + 10)