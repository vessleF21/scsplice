library(ggplot2)
library(dplyr)
library(tidyr)

df <- data.frame(
  Cis_eGene = c("PABPC1L", "AMPH", "WEE2-AS1", "MGAT5", "ZFP90", "ABLIM1", "IVNS1ABP"),
  Trans_sGene = c("NT5DC1", "MAN1A2", "SNRPG", "TWEM167A", "FUOM", "UCHL5", "RPS24"),
  Cell_Type = c("BNaive", "CD4EM", "CD4EM", "CD4Naive", "CD8GZMK", "NKDim", "NKDim"),
  H4 = c(0.87, 0.84, 0.89, 0.77, 0.99, 0.98, 1.00)
)

heat_data <- df %>%
  mutate(Interaction = paste(Cis_eGene, Trans_sGene, sep = "-")) %>%
  select(Interaction, Cell_Type, H4)

ggplot(heat_data, aes(x = Cell_Type, y = Interaction)) +
  geom_tile(aes(fill = H4), color = "white", size = 0.8) +
  geom_text(aes(label = sprintf("%.2f", H4)), color = "white", size = 4, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
    limits = c(0.7, 1.0),
    name = "H4",
    guide = guide_colorbar(
      barwidth = 15, 
      barheight = 0.8,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, color = "gray30", margin = margin(b = 10)),
    plot.caption = element_text(color = "gray40", hjust = 0.5, margin = margin(t = 10))
  ) +
  labs(
    x = "Cell_Type",
    y = "Cis-eGene - Trans-sGene",
    title = "eQTL-sQTL H4",
    subtitle = "Sig",
    caption = "eQTL/sQTL"
  ) +
  coord_fixed(ratio = 0.8)


 ggsave("eqtl_sqtl_heatmap.pdf", width = 12, height = 8, dpi = 300)
