

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(gridExtra)
})

CT_ORDER_F <- c(
  "BMem",     "CD4Naive", "CD4EM",
  "CD8Naive", "CD8GZMH",  "CD8GZMK",
  "MAIT",     "NKBright",  "cM", "ncM"
)

PATHWAY_ORDER <- c(
  "Retrograde endocannabinoid signaling",
  "Thermogenesis",
  "Spinocerebellar ataxia",
  "Alzheimer disease",
  "Prion disease",
  "Apoptosis",
  "Endocytosis",
  "Lysosome",
  "Natural killer cell mediated cytotoxicity",
  "B cell receptor signaling pathway"
)

.read_metascape <- function(ct, type) {
  fpath <- paste0("metascape/", type, "/", ct, "_enrichment.txt")
  if (!file.exists(fpath)) {
    return(data.frame(Description = PATHWAY_ORDER, neg_log10_p = 0))
  }
  df <- read.table(fpath, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE) |>
    dplyr::select(Description, log10p = `Log10(P)`) |>
    dplyr::mutate(neg_log10_p = -log10p) |>
    dplyr::filter(Description %in% PATHWAY_ORDER)
  missing <- setdiff(PATHWAY_ORDER, df$Description)
  if (length(missing) > 0) {
    df <- rbind(df, data.frame(Description = missing, log10p = 0, neg_log10_p = 0))
  }
  df
}

enrich_list <- list()
for (ct in CT_ORDER_F) {
  res_e <- .read_metascape(ct, "egene")
  res_s <- .read_metascape(ct, "sgene")
  enrich_list[[ct]] <- data.frame(
    cell_type         = ct,
    pathway           = PATHWAY_ORDER,
    neg_log10_p_egene = res_e$neg_log10_p[match(PATHWAY_ORDER, res_e$Description)],
    neg_log10_p_sgene = res_s$neg_log10_p[match(PATHWAY_ORDER, res_s$Description)]
  )
}

heatmap_raw <- do.call(rbind, enrich_list)
heatmap_raw[is.na(heatmap_raw)] <- 0

write.table(heatmap_raw, "panel_f_data.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

.build_mat <- function(dat, val_col) {
  mat <- dat |>
    dplyr::select(pathway, cell_type, !!val_col) |>
    tidyr::pivot_wider(names_from = cell_type, values_from = !!val_col) |>
    tibble::column_to_rownames("pathway") |>
    as.matrix()
  mat <- mat[PATHWAY_ORDER, CT_ORDER_F]
  mat[is.na(mat)] <- 0
  mat
}

mat_egene <- .build_mat(heatmap_raw, "neg_log10_p_egene")
mat_sgene <- .build_mat(heatmap_raw, "neg_log10_p_sgene")

.draw_heatmap <- function(mat, colors, breaks, title) {
  pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color        = colors,
    breaks       = breaks,
    cellwidth    = 14,
    cellheight   = 14,
    border_color = "grey90",
    fontsize_row = 8,
    fontsize_col = 8,
    main         = title,
    silent       = TRUE
  )
}

fig_f_egene <- .draw_heatmap(
  mat    = mat_egene,
  colors = colorRampPalette(c("white", "#2166AC"))(50),
  breaks = seq(0, 5,  length.out = 51),
  title  = "Cis-eGene-related pathway -log\u2081\u2080(P)"
)

fig_f_sgene <- .draw_heatmap(
  mat    = mat_sgene,
  colors = colorRampPalette(c("white", "#B2182B"))(50),
  breaks = seq(0, 18, length.out = 51),
  title  = "Cis-sGene-related pathway -log\u2081\u2080(P)"
)

pdf("panel_f.pdf", width = 16, height = 6)
gridExtra::grid.arrange(fig_f_egene$gtable, fig_f_sgene$gtable, ncol = 2)
dev.off()


