

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(fs)
  library(stringr)
  library(dplyr)
  library(ggpmisc)
})

OUTPUT_SCATTER  <- "square_plot1.pdf"
OUTPUT_TABLE    <- "multi_sqtl_gene.txt"
POINT_SIZE      <- 4
POINT_STROKE    <- 2
LINE_WIDTH      <- 1.5
ANNO_SIZE       <- 5
AXIS_TEXT_SIZE  <- 30
AXIS_TITLE_SIZE <- 30

ct_color_pal <- c(
  "BNav"     = "#01579B",
  "BMem"     = "#A5EDFF",
  "CD4EM"    = "#1B5E20",
  "CD4Naive" = "#8BC34A",
  "CD4Treg"  = "#D4E157",
  "CD8GZMH"  = "#FF5722",
  "CD8GZMK"  = "#FFCDD2",
  "CD8Naive" = "#F9A825",
  "MAIT"     = "#FFEE58",
  "MonocM"   = "#6A1B9A",
  "NoClaM"   = "#E1BEE7",
  "NKBright" = "#82581F",
  "NKDim"    = "#D7CCC8"
)

ct_color_table <- tibble(
  ct_name = c(
    "BMem",
    "BNav",
    "CD4EM",
    "CD4Naive",
    "CD4Treg",
    "CD8GZMH",
    "CD8GZMK",
    "CD8Naive",
    "MAIT",
    "MonocM",
    "NKBright",
    "NKDim",
    "NoClaM"
  ),
  hex_color = c(
    "#ffffcc",
    "#ffeda0",
    "#bf812d",
    "#efedf5",
    "#dadaeb",
    "#bcbddc",
    "#9e9ac8",
    "#807dba",
    "#8c510a",
    "#74a9cf",
    "#3690c0",
    "#a6bddb",
    "#74c476"
  )
)

raw_input <- read.table(
  "data1.txt",
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE,
  comment.char     = ""
)

ct_multi_sqtl <- raw_input |>
  dplyr::group_by(cell_type) |>
  dplyr::summarise(
    n_multi_sqtl = sum(sqtl > 1),
    n_total      = dplyr::n(),
    prop_multi   = n_multi_sqtl / n_total,
    .groups      = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(prop_multi)) |>
  dplyr::left_join(ct_color_table, by = c("cell_type" = "ct_name"))

print(ct_multi_sqtl)

write.table(ct_multi_sqtl, file = OUTPUT_TABLE,
            row.names = FALSE, quote = FALSE, sep = "\t")

ct_multi_sqtl <- read.table(
  OUTPUT_TABLE,
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE,
  comment.char     = ""
)

pearson_res <- cor.test(
  ct_multi_sqtl$n_total,
  ct_multi_sqtl$prop_multi,
  method = "pearson"
)

corr_label <- paste0(
  "Pearson's r = ", round(pearson_res$estimate, 2),
  "\nP = ", formatC(pearson_res$p.value, format = "e", digits = 2)
)

fig_scatter <- ggplot(ct_multi_sqtl, aes(x = n_total, y = prop_multi)) +
  geom_point(aes(color = cell_type), size = POINT_SIZE,
             shape = 21, fill = "white", stroke = POINT_STROKE) +
  scale_color_manual(values = ct_color_pal) +
  stat_poly_eq(
    formula  = y ~ x,
    aes(label = after_stat(eq.label)),
    size = 4, label.x = 3, label.y = 0.1, color = "darkred"
  ) +
  stat_poly_line(formula = y ~ x, color = "darkblue", linewidth = LINE_WIDTH) +
  ggplot2::annotate(
    geom  = "text",
    x = Inf, y = Inf, label = corr_label,
    hjust = 1.1, vjust = 1.5, size = ANNO_SIZE, color = "darkgreen"
  ) +
  labs(
    x = "Number of sGenes",
    y = "Proportion of sGenes\nwith >1 sSNPs"
  ) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title       = element_text(face = "bold", size = AXIS_TITLE_SIZE, color = "black"),
    axis.text        = element_text(size = AXIS_TEXT_SIZE, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    axis.line        = element_line(color = "black", linewidth = 0.6)
  )

print(fig_scatter)
ggsave(OUTPUT_SCATTER, fig_scatter, width = 7, height = 7)


