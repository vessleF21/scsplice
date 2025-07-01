rm(list = ls())
gc()

library(dplyr)
library(ggplot2)
library(ggpubr)

dat.sig <- data.frame(
  celltype = c("BNav", "BNav", "BNav", "CD4EM", "CD4Naive", "CD4Naive", "CD4Treg", "CD4Treg", "CD4Treg", "CD8GZMH", "CD8GZMH", "CD8GZMK", "CD8GZMK", "CD8GZMK", "CD8Naive", "CD8Naive", "CD8Naive", "MAIT", "MonocM", "MonocM", "NKDim", "NKDim", "NoClaM", "NoClaM", "NoClaM"),
  country = c("age_45_65vs45low", "age_high65vs45_65", "age_high65vs45low", "age_45_65vs45low",
              "age_45_65vs45low", "age_high65vs45_65", "age_45_65vs45low", "age_high65vs45_65",
              "age_high65vs45low", "age_high65vs45_65", "age_high65vs45low", "age_45_65vs45low",
              "age_high65vs45_65", "age_high65vs45low", "age_45_65vs45low", "age_high65vs45_65",
              "age_high65vs45low", "age_high65vs45_65", "age_high65vs45_65", "age_high65vs45low",
              "age_high65vs45_65", "age_high65vs45low", "age_45_65vs45low", "age_high65vs45_65",
              "age_high65vs45low"),
  count = c(1, 4, 5, 4, 3, 7, 3, 7, 3, 2, 5, 1, 1, 7, 9, 6, 33, 2, 2, 3, 9, 21, 1, 2, 3)
)

color_mapping <- c(
  "BNav" = "#01579B",
  "BMem" = "#A5EDFF",
  "CD4EM" = "#1B5E20",
  "CD4Naive" = "#8BC34A",
  "CD4Treg" = "#D4E157",
  "CD8GZMH" = "#FF5722",
  "CD8GZMK" = "#FFCDD2",
  "CD8Naive" = "#F9A825",
  "MAIT" = "#FFEE58",
  "MonocM" = "#6A1B9A",
  "NoClaM" = "#E1BEE7",
  "NKBright" = "#82581F",
  "NKDim" = "#D7CCC8"
)

dat.sig.sum <- dat.sig %>%
  group_by(celltype) %>%
  summarise(total_count = sum(count), .groups = 'drop')

celltype_order <- c("BNav", "CD4Naive", "CD4EM", "CD4Treg", "CD8Naive", "CD8GZMH", "CD8GZMK", "MAIT", "NKDim", "MonocM", "NoClaM")
dat.sig.sum$celltype <- factor(dat.sig.sum$celltype, levels = celltype_order)
dat.sig$celltype <- factor(dat.sig$celltype, levels = celltype_order)

dat.sig.percent <- dat.sig %>%
  group_by(celltype) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(percent = count / total * 100)

country_pal <- c("#f06152", "#76afda", "#7fb961")
names(country_pal) <- c('age_45_65vs45low', 'age_high65vs45_65', 'age_high65vs45low')

f2 <- ggplot(dat.sig.sum, aes(x = celltype, y = total_count, fill = celltype)) +
  geom_bar(stat = 'identity', color = 'gray', width = 0.8) +
  ylab('Total Count') + xlab('Cell Type') +
  geom_text(aes(label = total_count), vjust = -0.3, size = 4, color = 'black', fontface = 'bold') +
  theme_minimal() +
  scale_fill_manual(values = color_mapping) +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 12, face = "bold", color = '#404040'),
    axis.text.y = element_text(size = 12, face = "bold", color = '#404040')
  )

f3 <- ggplot(dat.sig.percent, aes(x = celltype, y = percent, fill = country)) +
  geom_bar(stat = "identity", position = "stack", color = "white", width = 0.8) +
  geom_text(aes(label = paste0(round(percent), "%")), position = position_stack(vjust = 0.5), size = 4, color = "black") +
  scale_fill_manual(values = country_pal) +
  labs(x = "Cell Type", y = "") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", size = 1),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 12, face = "bold", color = '#404040'),
    axis.text.y = element_text(size = 12, face = "bold", color = '#404040')
  )

f23 <- ggarrange(f2, f3, nrow = 2, align = 'v', heights = c(2, 3), common.legend = TRUE)

ggsave("f2_plot.pdf", plot = f2, width = 10, height = 6)
ggsave("f3_plot.pdf", plot = f3, width = 10, height = 6)
ggsave('Fig.23f.pdf', plot = f23, width = 10, height = 12)
