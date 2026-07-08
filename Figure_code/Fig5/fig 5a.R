


rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

CT_ORDER <- c(
  "BNav", "BMem", "CD4Naive", "CD4EM", "CD4Treg",
  "CD8Naive", "CD8GZMH", "CD8GZMK", "MAIT",
  "NKDim", "NKBright", "cM", "ncM"
)
COL_BLUE   <- "#4472C4"
COL_ORANGE <- "#ED7D31"

bar_list <- list()
for (ct in CT_ORDER) {
  egene <- read.table(paste0("egene/", ct, "_egene.txt"),
                      header = TRUE, stringsAsFactors = FALSE)
  sgene <- read.table(paste0("sgene/", ct, "_sgene.txt"),
                      header = TRUE, stringsAsFactors = FALSE)
  n_total   <- nrow(egene)
  n_overlap <- length(intersect(egene$gene, sgene$gene))
  bar_list[[ct]] <- data.frame(
    cell_type = ct,
    n_total   = n_total,
    n_overlap = n_overlap,
    n_only    = n_total - n_overlap
  )
}

bar_raw <- do.call(rbind, bar_list)
write.table(bar_raw, "panel_a_data.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

bar_data <- bar_raw |>
  dplyr::mutate(
    pct_overlap = n_overlap / n_total * 100,
    pct_only    = n_only    / n_total * 100,
    label       = paste0(cell_type, "\n(n=", n_total, ")")
  ) |>
  tidyr::pivot_longer(
    cols      = c(pct_overlap, pct_only),
    names_to  = "category",
    values_to = "percentage"
  ) |>
  dplyr::mutate(
    category  = factor(category,
                       levels = c("pct_overlap", "pct_only"),
                       labels = c("Cis-eGene-Cis-sGene overlap", "Cis-eGene only")),
    cell_type = factor(cell_type, levels = CT_ORDER)
  )

fig_a <- ggplot(bar_data, aes(x = cell_type, y = percentage, fill = category)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_manual(values = c(COL_BLUE, COL_ORANGE)) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = paste0(seq(0, 100, 25), "%")
  ) +
  scale_x_discrete(labels = setNames(bar_data$label, bar_data$cell_type)) +
  labs(x = NULL, y = "Percentage (%)", fill = NULL) +
  theme_classic() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y     = element_text(size = 8),
    axis.title.y    = element_text(size = 9),
    legend.position = "top",
    legend.text     = element_text(size = 8),
    panel.grid      = element_blank()
  )

ggsave("panel_a.pdf", fig_a, width = 10, height = 5, dpi = 300)



