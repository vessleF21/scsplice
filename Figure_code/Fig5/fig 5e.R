

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

COL_TEAL     <- "#70AD47"
COL_RED      <- "#FF0000"
LD_THRESHOLD <- 0.8

partial_genes <- read.table(
  "panel_d_data.txt",
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE
) |>
  dplyr::filter(category == "Partial overlap") |>
  dplyr::pull(gene)

r2_list <- lapply(partial_genes, function(g) {
  egene_snp <- read.table(paste0("egene_snp/", g, ".txt"),
                          header = TRUE, stringsAsFactors = FALSE)$lead_snp[1]
  sgene_snp <- read.table(paste0("sgene_snp/", g, ".txt"),
                          header = TRUE, stringsAsFactors = FALSE)$lead_snp[1]
  ld_file   <- read.table(paste0("ld/", g, ".ld"),
                          header = TRUE, stringsAsFactors = FALSE)
  r2_val <- ld_file |>
    dplyr::filter(
      (SNP_A == egene_snp & SNP_B == sgene_snp) |
      (SNP_A == sgene_snp & SNP_B == egene_snp)
    ) |>
    dplyr::pull(R2)
  data.frame(
    gene     = g,
    r2       = ifelse(length(r2_val) > 0, r2_val[1], 0),
    category = ifelse(
      length(r2_val) > 0 && r2_val[1] >= LD_THRESHOLD,
      "R\u00b2 \u2265 0.8",
      "R\u00b2 < 0.8"
    )
  )
})

pie_e_raw <- do.call(rbind, r2_list)
write.table(pie_e_raw, "panel_e_data.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

pie_e_df <- pie_e_raw |>
  dplyr::count(category) |>
  dplyr::mutate(
    pct      = n / sum(n) * 100,
    label    = paste0(n, " (", round(pct), "%)"),
    category = factor(category, levels = c("R\u00b2 < 0.8", "R\u00b2 \u2265 0.8"))
  )

fig_e <- ggplot(pie_e_df, aes(x = "", y = n, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(COL_RED, COL_TEAL), name = NULL) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  theme_void() +
  theme(legend.text = element_text(size = 8))

ggsave("panel_e.pdf", fig_e, width = 5, height = 5, dpi = 300)

