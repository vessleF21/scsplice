

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(VennDiagram)
  library(grid)
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

all_sgene <- unique(unlist(lapply(CT_ORDER, function(ct) {
  read.table(paste0("sgene/", ct, "_sgene.txt"),
             header = TRUE, stringsAsFactors = FALSE)$gene
})))

N_EGENE  <- length(all_egene)
N_SGENE  <- length(all_sgene)
N_SHARED <- length(intersect(all_egene, all_sgene))

venn_summary <- data.frame(
  category = c("cis-eGene", "cis-sGene", "shared"),
  count    = c(N_EGENE, N_SGENE, N_SHARED)
)
write.table(venn_summary, "panel_b_data.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

fig_b <- VennDiagram::draw.pairwise.venn(
  area1        = N_EGENE,
  area2        = N_SGENE,
  cross.area   = N_SHARED,
  category     = c("Cis-eGene", "Cis-sGene"),
  fill         = c("grey70", "grey40"),
  alpha        = 0.5,
  cex          = 1.2,
  cat.cex      = 1.0,
  cat.fontface = "italic",
  euler.d      = TRUE,
  scaled       = TRUE
)

pdf("panel_b.pdf", width = 5, height = 5)
grid::grid.draw(fig_b)
dev.off()

