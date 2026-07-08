


rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

COL_TEAL <- "#70AD47"
COL_RED  <- "#FF0000"
COL_PINK <- "#FFC0CB"

shared_genes <- readLines("shared_genes.txt")

category_list <- lapply(shared_genes, function(g) {
  egene_snps <- read.table(paste0("egene_snp/", g, ".txt"),
                           header = TRUE, stringsAsFactors = FALSE)$snp
  sgene_snps <- read.table(paste0("sgene_snp/", g, ".txt"),
                           header = TRUE, stringsAsFactors = FALSE)$snp
  if (length(intersect(egene_snps, sgene_snps)) == 0) {
    "Completely distinct sets"
  } else if (setequal(egene_snps, sgene_snps)) {
    "Identical sets"
  } else {
    "Partial overlap"
  }
})

pie_d_raw <- data.frame(gene = shared_genes, category = unlist(category_list))
write.table(pie_d_raw, "panel_d_data.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

pie_d_df <- pie_d_raw |>
  dplyr::count(category) |>
  dplyr::mutate(
    pct      = n / sum(n) * 100,
    label    = paste0(n, " (", round(pct, 2), "%)"),
    category = factor(category,
                      levels = c("Completely distinct sets",
                                 "Partial overlap",
                                 "Identical sets"))
  )

fig_d <- ggplot(pie_d_df, aes(x = "", y = n, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(COL_TEAL, COL_RED, COL_PINK), name = NULL) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  theme_void() +
  theme(legend.text = element_text(size = 8))

ggsave("panel_d.pdf", fig_d, width = 5, height = 5, dpi = 300)



