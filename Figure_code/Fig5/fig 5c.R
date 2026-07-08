



rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

CT_ORDER <- c(
  "BNav", "BMem", "CD4Naive", "CD4EM", "CD4Treg",
  "CD8Naive", "CD8GZMH", "CD8GZMK", "MAIT",
  "NKDim", "NKBright", "cM", "ncM"
)

all_egene <- unique(unlist(lapply(CT_ORDER, function(ct) {
  read.table(paste0("egene/", ct, "_egene.txt"),
             header = TRUE, stringsAsFactors = FALSE)$gene
})))

all_sgene_per_ct <- lapply(CT_ORDER, function(ct) {
  read.table(paste0("sgene/", ct, "_sgene.txt"),
             header = TRUE, stringsAsFactors = FALSE)$gene
})
names(all_sgene_per_ct) <- CT_ORDER

shared_genes <- intersect(all_egene, unique(unlist(all_sgene_per_ct)))

n_ct_per_gene <- sapply(shared_genes, function(g) {
  sum(sapply(all_sgene_per_ct, function(sg) g %in% sg))
})

pie_c_raw <- as.data.frame(table(n_ct_per_gene)) |>
  dplyr::rename(n_celltypes = n_ct_per_gene, n_genes = Freq) |>
  dplyr::mutate(n_celltypes = as.integer(as.character(n_celltypes)))

write.table(pie_c_raw, "panel_c_data.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

pie_c_data <- pie_c_raw |>
  dplyr::arrange(n_celltypes) |>
  dplyr::mutate(
    pct   = n_genes / sum(n_genes) * 100,
    label = ifelse(pct > 3, paste0(n_genes), "")
  )

fig_c <- ggplot(pie_c_data,
                aes(x = "", y = n_genes, fill = factor(n_celltypes))) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  scale_fill_manual(
    values = colorRampPalette(c("#E8F5E9", "#1B5E20"))(13),
    name   = "Number of shared\ncell types"
  ) +
  annotate("text", x = 0, y = 0,
           label    = paste0(length(shared_genes), "\noverlapped\ngenes"),
           size     = 3.5,
           fontface = "bold") +
  theme_void() +
  theme(legend.text = element_text(size = 8))

ggsave("panel_c.pdf", fig_c, width = 5, height = 5, dpi = 300)

