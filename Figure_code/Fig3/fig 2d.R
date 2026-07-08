


library(ggplot2)
library(dplyr)
library(patchwork)

raw <- read.delim("onek1k.bed", header = FALSE)

gene_cov <- raw |>
  dplyr::group_by(V5) |>
  dplyr::summarise(
    all_bases     = sum(V9),
    covered_bases = sum(V8),
    sumReads      = sum(V7),
    .groups       = "drop"
  ) |>
  dplyr::mutate(perc_covered = covered_bases / all_bases)

gene_cov <- gene_cov |>
  dplyr::filter(sumReads > 5) |>
  na.omit()

bin_levels <- c("1 - 1000", "1000 - 10000", "10000 - 100000", "100000+")

gene_cov <- gene_cov |>
  dplyr::mutate(
    bin = dplyr::case_when(
      sumReads <= 1000   ~ "1 - 1000",
      sumReads <= 10000  ~ "1000 - 10000",
      sumReads <= 100000 ~ "10000 - 100000",
      TRUE               ~ "100000+"
    ),
    bin = factor(bin, levels = bin_levels)
  )

pal_gray <- c(
  "1 - 1000"       = "#FFFFFF",
  "1000 - 10000"   = "#C8C8C8",
  "10000 - 100000" = "#888888",
  "100000+"        = "#444444"
)

med <- median(gene_cov$perc_covered)
cat(sprintf("Overall median: %.3f\n", med))

fig_left <-
  ggplot(gene_cov, aes(x = bin, y = perc_covered, fill = bin)) +
  geom_boxplot(
    width         = 0.65,
    color         = "black",
    linewidth     = 0.5,
    outlier.shape = NA
  ) +
  scale_fill_manual(values = pal_gray, name = "Read Number") +
  scale_y_continuous(
    limits = c(0, 1.02),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0", "0.25", "0.50", "0.75", "1.00")
  ) +
  labs(y = "Fraction of\nbase coverage", x = "Number of reads per gene") +
  theme_classic() +
  theme(
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.title.y      = element_text(size = 10),
    legend.position   = c(0.38, 0.35),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size   = unit(0.45, "cm"),
    legend.key        = element_rect(color = "black", linewidth = 0.4),
    legend.text       = element_text(size = 8),
    legend.title      = element_text(size = 9)
  )

gene_cov_all       <- gene_cov
gene_cov_all$group <- "All"

fig_right <-
  ggplot(gene_cov_all, aes(x = group, y = perc_covered)) +
  geom_violin(fill = "#D5D0BC", color = "black", trim = FALSE, linewidth = 0.5) +
  geom_boxplot(
    width         = 0.12,
    fill          = "#8B8B6B",
    color         = "black",
    linewidth     = 0.5,
    outlier.shape = NA
  ) +
  geom_hline(yintercept = med, color = "red", linetype = "dashed", linewidth = 0.7) +
  ggplot2::annotate(
    geom  = "text",
    x     = 1.35,
    y     = med + 0.04,
    label = sprintf("%.3f", med),
    color = "red", size = 3.5, hjust = 0
  ) +
  scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
  theme_classic() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title   = element_blank(),
    axis.text.x  = element_text(size = 9)
  )

fig_final <- fig_left + fig_right + plot_layout(widths = c(3, 1))

print(fig_final)

ggsave("coverage_plot.pdf", fig_final, width = 10, height = 5, dpi = 300, device = cairo_pdf)
ggsave("coverage_plot.png", fig_final, width = 10, height = 5, dpi = 300)