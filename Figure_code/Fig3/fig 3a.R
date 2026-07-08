

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(fs)
  library(stringr)
  library(dplyr)
  library(ggpmisc)
})

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

ref_gene <- read_gene_annotation(gene_annotation_path)

sample_table <- read_donor_list(
  donor_list_dir,
  cell_types = cell_types
)

ct_stats <- sample_table |>
  dplyr::group_by(cell_type) |>
  dplyr::summarise(
    n_donors     = dplyr::n(),
    mean_cells   = mean(number),
    median_cells = median(number),
    total_cells  = sum(number),
    .groups      = "drop"
  )

print(ct_stats)

raw_perm <- read_perm_qtl_all_cell_types(cell_types, data_dir = perm_dir)
cat(sprintf("raw_perm dim: %d x %d\n", nrow(raw_perm), ncol(raw_perm)))

expanded_perm <- split_multiple_genes(raw_perm)
cat(sprintf("expanded_perm dim: %d x %d\n", nrow(expanded_perm), ncol(expanded_perm)))

annot_perm <- expanded_perm |>
  dplyr::left_join(ref_gene[, c("gene_name", "gene_id", "type")], by = "gene_name")
cat(sprintf("annot_perm dim: %d x %d\n", nrow(annot_perm), ncol(annot_perm)))

raw_cond <- read_indep_qtl_all_cell_types(cell_types, conditional_dir)
cat(sprintf("raw_cond dim: %d x %d\n", nrow(raw_cond), ncol(raw_cond)))

expanded_cond <- split_multiple_genes(raw_cond)
cat(sprintf("expanded_cond dim: %d x %d\n", nrow(expanded_cond), ncol(expanded_cond)))

annot_cond <- expanded_cond |>
  dplyr::left_join(ref_gene[, c("gene_name", "gene_id", "type")], by = "gene_name")
cat(sprintf("annot_cond dim: %d x %d\n", nrow(annot_cond), ncol(annot_cond)))

n_unique <- annot_cond |>
  dplyr::select(gene_name, variant_id, cell_type) |>
  dplyr::distinct() |>
  nrow()
cat(sprintf("Unique (gene, variant, celltype): %d\n", n_unique))

dedup_cond <- annot_cond |>
  dplyr::select(-gene_id) |>
  dplyr::distinct()
cat(sprintf("dedup_cond after dedup: %d rows\n", nrow(dedup_cond)))

annot_cond |>
  dplyr::select(gene_name, type) |>
  dplyr::distinct() |>
  dplyr::count(type, sort = TRUE) |>
  print()

annot_perm |>
  dplyr::select(gene_name, type) |>
  dplyr::distinct() |>
  dplyr::count(type, sort = TRUE) |>
  print()

gene_sqtl_summary <- annot_cond |>
  dplyr::group_by(cell_type, gene_name) |>
  dplyr::summarise(sqtl = dplyr::n(), .groups = "drop")

write.table(gene_sqtl_summary, file = "data1.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

cat(sprintf("rows: %d, cols: %d\n", nrow(gene_sqtl_summary), ncol(gene_sqtl_summary)))

DIAMOND_COLOR <- "#6A3D9A"
DIAMOND_RANGE <- c(3, 12)
BAR_WIDTH     <- 0.85
BAR_LINEWIDTH <- 0.3
LABEL_SIZE    <- 3.2
LABEL_VJUST   <- -1.8
BASE_FONT     <- 12
OUTPUT_DIR    <- "sqtl_summary12/"

fill_pal <- c(
  "1"   = "#FDE725",
  "2"   = "#7AD151",
  ">=3" = "#22A884"
)

raw_input <- read.table(
  "data1.txt",
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE,
  comment.char     = ""
)

binned_sqtl <- raw_input |>
  dplyr::mutate(
    sqtl = dplyr::if_else(sqtl >= 3, ">=3", as.character(sqtl)),
    sqtl = factor(sqtl, levels = c("1", "2", ">=3"))
  )

ct_summary <- binned_sqtl |>
  dplyr::group_by(cell_type) |>
  dplyr::summarise(n_sgene = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(n_sgene))

ct_order   <- dplyr::pull(ct_summary, cell_type)
norm_scale <- max(ct_summary$n_sgene) * 1.10

dir.create(OUTPUT_DIR, showWarnings = FALSE)

fig_sqtl <- ggplot() +
  geom_bar(
    data      = binned_sqtl,
    aes(x = cell_type, fill = sqtl),
    position  = "fill",
    color     = "black",
    linewidth = BAR_LINEWIDTH,
    width     = BAR_WIDTH
  ) +
  geom_point(
    data  = ct_summary,
    aes(x = cell_type, y = n_sgene / norm_scale, size = n_sgene),
    color = DIAMOND_COLOR,
    shape = 18
  ) +
  geom_text(
    data  = ct_summary,
    aes(x = cell_type, y = n_sgene / norm_scale, label = n_sgene),
    color    = DIAMOND_COLOR,
    vjust    = LABEL_VJUST,
    size     = LABEL_SIZE,
    fontface = "bold"
  ) +
  scale_y_continuous(
    name     = "Proportion of sGenes",
    expand   = c(0, 0),
    labels   = scales::percent_format(),
    sec.axis = sec_axis(
      trans  = ~ . * norm_scale,
      name   = "Number of sGenes",
      breaks = scales::pretty_breaks(n = 6)
    )
  ) +
  scale_x_discrete(limits = ct_order) +
  scale_fill_manual(
    values = fill_pal,
    name   = "sQTL/sGene",
    labels = c("1", "2", "\u22653")
  ) +
  scale_size_continuous(
    name   = "Number of sGenes",
    range  = DIAMOND_RANGE,
    breaks = scales::pretty_breaks(n = 5)
  ) +
  guides(
    fill = guide_legend(order = 1, title = "sQTL/sGene"),
    size = guide_legend(
      order        = 2,
      title        = "Number of sGenes",
      override.aes = list(shape = 18, color = DIAMOND_COLOR)
    )
  ) +
  theme_minimal(base_size = BASE_FONT) +
  theme(
    axis.title         = element_text(face = "bold"),
    axis.text.x        = element_text(angle = 60, hjust = 1,
                                      vjust = 1, color = "black", size = 10),
    axis.text.y        = element_text(color = "black", size = 10),
    axis.text.y.right  = element_text(color = DIAMOND_COLOR, size = 10),
    axis.title.y.right = element_text(color = DIAMOND_COLOR, margin = margin(l = 10)),
    axis.line          = element_line(color = "black"),
    axis.line.y.right  = element_line(color = DIAMOND_COLOR),
    axis.ticks.y.right = element_line(color = DIAMOND_COLOR),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position    = "bottom",
    legend.title       = element_text(face = "bold"),
    legend.text        = element_text(size = 9),
    legend.key.size    = unit(0.8, "cm")
  )

print(fig_sqtl)

ggsave(file.path(OUTPUT_DIR, "sqtl_barplot.pdf"), fig_sqtl, width = 16, height = 10)
ggsave(file.path(OUTPUT_DIR, "sqtl_barplot.png"), fig_sqtl, width = 16, height = 10, dpi = 300)



