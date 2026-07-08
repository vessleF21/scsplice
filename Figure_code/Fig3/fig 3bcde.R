


rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(data.table)
  library(dplyr)
  library(ggpmisc)
})

donor_list_dir       <- "./split1/"
gene_annotation_path <- "gencode.v32.primary_assembly.annotation.gtf.genelist.tsv"
conditional_dir      <- "chr_top/"
numerator_file       <- "merged_oneK1K_perind_numers.counts.gz"

.OUT_SUMMARY <- "cell_type_summary_1.txt"
.OUT_DONOR   <- "donor_list_2.txt"
.OUT_DIR     <- "sqtl_summary/"

dir.create(.OUT_DIR, showWarnings = FALSE)

CT_PALETTE <- c(
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

CT_COLOR_REF <- tibble(
  cell_type = names(CT_PALETTE),
  hex_color = unname(CT_PALETTE)
)

gene_ref <- read_gene_annotation(gene_annotation_path)

cond_raw      <- read_indep_qtl_all_cell_types(cell_types, conditional_dir)
cond_expanded <- split_multiple_genes(cond_raw)

cond_annot <- cond_expanded |>
  dplyr::left_join(
    gene_ref[, c("gene_name", "gene_id", "type")],
    by = "gene_name"
  )

sqtl_gene <- cond_annot |>
  dplyr::group_by(cell_type, gene_name) |>
  dplyr::summarise(
    n_sqtl    = dplyr::n(),
    max_slope = max(abs(fore_slope)),
    .groups   = "drop"
  )

sgene_ct <- sqtl_gene |>
  dplyr::group_by(cell_type) |>
  dplyr::summarise(
    n_sgene      = dplyr::n(),
    mean_slope   = mean(max_slope),
    median_slope = median(max_slope),
    .groups      = "drop"
  )

donor_raw <- read_donor_list(donor_list_dir, cell_types = cell_types)

numer_raw <- fread(numerator_file)
colnames(numer_raw)[1] <- "chrom"

col_info <- data.table(col_names = colnames(numer_raw)[-1]) |>
  dplyr::mutate(
    DCP_ID    = sapply(strsplit(col_names, "@"), \(x) x[length(x) - 1]),
    cell_type = sapply(
      strsplit(
        sapply(strsplit(col_names, "@"), \(x) x[length(x)]),
        "\\."
      ),
      \(x) x[1]
    )
  )

lib_raw <- data.frame(colSums(numer_raw[, -1])) |>
  as.data.table(keep.rownames = TRUE)
colnames(lib_raw) <- c("col_names", "lib_size")

lib_annot <- lib_raw |>
  dplyr::left_join(col_info, by = "col_names")

donor_lib <- donor_raw |>
  dplyr::left_join(lib_annot, by = c("DCP_ID", "cell_type"))

ct_base <- donor_lib |>
  dplyr::group_by(cell_type) |>
  dplyr::summarise(
    n_donors       = dplyr::n(),
    mean_cells     = mean(number),
    median_cells   = median(number),
    total_cells    = sum(number),
    mean_libsize   = mean(lib_size),
    median_libsize = median(lib_size),
    total_libsize  = sum(lib_size),
    .groups        = "drop"
  )

ct_final <- ct_base |>
  dplyr::left_join(sgene_ct,     by = "cell_type") |>
  dplyr::left_join(CT_COLOR_REF, by = "cell_type") |>
  dplyr::mutate(color = CT_PALETTE[cell_type])

.save_tsv <- function(x, path) {
  write.table(x, file = path, row.names = FALSE, quote = FALSE, sep = "\t")
  cat(sprintf("saved: %s\n", path))
}

.save_tsv(ct_final,  .OUT_SUMMARY)
.save_tsv(donor_lib, .OUT_DONOR)

.base_scatter <- function(dat, x_var, y_var,
                           x_lab, y_lab,
                           formula   = y ~ x,
                           log_y     = FALSE,
                           text_size = 14) {
  r_label <- sprintf("Pearson's r = %.2f",
                     cor(dat[[x_var]], dat[[y_var]], method = "pearson"))

  gg <- ggplot(dat, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(aes(color = cell_type), size = 5, shape = 21,
               fill = "white", stroke = 2) +
    scale_color_manual(values = CT_PALETTE) +
    stat_poly_eq(formula  = formula,
                 aes(label = after_stat(eq.label)),
                 size = 4, label.x = "right", color = "darkred") +
    stat_poly_line(formula = formula, color = "darkblue", linewidth = 1.5) +
    ggplot2::annotate("text", x = Inf, y = Inf, label = r_label,
                      hjust = 1.1, vjust = 1.5, size = 5, color = "darkgreen") +
    labs(x = x_lab, y = y_lab) +
    theme_classic() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.title       = element_text(face = "bold", size = 16, color = "black"),
      axis.text        = element_text(size = text_size, color = "black"),
      panel.grid       = element_blank(),
      legend.position  = "none",
      axis.line        = element_line(color = "black", linewidth = 0.6)
    )

  if (log_y) gg <- gg + scale_y_log10()
  gg
}

.save_fig <- function(fig, path, w = 7, h = 7) {
  print(fig)
  ggsave(path, fig, width = w, height = h)
  cat(sprintf("saved: %s\n", path))
}

.save_fig(.base_scatter(ct_final, "n_donors",    "n_sgene",
                         "Number of donors", "Number of sGenes",
                         log_y = TRUE, text_size = 22),
          "sgene_vs_donor.pdf")

.save_fig(.base_scatter(ct_final, "mean_cells",  "n_sgene",
                         "Cells per donor", "Number of sGenes",
                         text_size = 22),
          "sgene_vs_mean_cells.pdf")

.save_fig(.base_scatter(ct_final, "n_donors",    "mean_slope",
                         "Number of donors", "Mean absolute effect size",
                         formula = y ~ x + poly(x, 2)),
          "effect_size_vs_sample_size.pdf")

.save_fig(.base_scatter(ct_final, "mean_libsize", "mean_slope",
                         "Library size", "Mean absolute effect size"),
          "effect_size_vs_library_size.pdf")