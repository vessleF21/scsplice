

library(data.table)
library(ggplot2)
library(patchwork)

CELL_ORDER <- c(
  "BNav", "BMem", "CD4EM", "CD4Naive", "CD4Treg",
  "CD8Naive", "CD8GZMH", "CD8GZMK", "MAIT",
  "NKDim", "NKBright", "MonocM", "NoClaM"
)

COLORS_CELL <- c(
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

args   <- commandArgs(trailingOnly = TRUE)
SNP    <- args[1]
intron <- args[2]
gene   <- args[3]

genocol  <- read.table("geno.txt", sep = "\t", check.names = FALSE)
genotype <- read.table(paste0(gene, "/snp/", SNP, ".txt"),
                       sep = "\t", check.names = FALSE,
                       stringsAsFactors = FALSE, colClasses = "character")
colnames(genotype) <- genocol

REF <- genotype$REF[1]
ALT <- genotype$ALT[1]

genotype[10:ncol(genotype)] <- lapply(genotype[10:ncol(genotype)], function(col) {
  sapply(col, function(gt) {
    if (gt %in% c("0/1", "1/0")) paste0(REF, ALT)
    else if (gt == "0/0")        paste0(REF, REF)
    else if (gt == "1/1")        paste0(ALT, ALT)
    else NA
  })
})

plots <- list()

for (ct in CELL_ORDER) {
  color     <- COLORS_CELL[[ct]]
  file_path <- paste0(gene, "/phe/cellty/", gene, "/",
                      intron, "@", SNP, "/",
                      ct, "@", intron, "@", SNP, ".txt")

  if (!file.exists(file_path)) {
    plots[[length(plots) + 1]] <- ggplot() +
      labs(title = paste("No Data:", ct), x = "Genotype", y = "Phenotype") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    next
  }

  col_phe   <- as.character(unlist(read.table(
    paste0(ct, "_phename.txt"), sep = " ", check.names = FALSE)))
  phenotype <- as.data.frame(fread(file_path, sep = " "))
  colnames(phenotype) <- col_phe
  genotype  <- as.data.frame(genotype)

  common_samples <- intersect(colnames(genotype), colnames(phenotype))[-1]

  merged <- as.data.frame(t(rbind(
    phenotype[, common_samples, drop = FALSE],
    genotype[,  common_samples, drop = FALSE]
  )))
  rownames(merged) <- NULL
  colnames(merged) <- c("phe", "geno")
  merged$geno <- factor(merged$geno,
                        levels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)))
  merged$phe  <- as.numeric(merged$phe)

  plots[[length(plots) + 1]] <- ggplot(merged, aes(x = geno, y = phe)) +
    geom_violin(trim = FALSE, alpha = 0, scale = "width",
                linewidth = 0.5, color = color) +
    geom_boxplot(width = 0.1, fill = "white",
                 outlier.shape = NA, color = color) +
    labs(title = ct, x = "Genotype", y = "Phenotype") +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0.5),
      legend.position   = "none",
      axis.line.x       = element_line(color = "black", linewidth = 0.5),
      axis.line.y       = element_line(color = "black", linewidth = 0.5),
      axis.ticks        = element_line(color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      panel.grid        = element_blank()
    )
}

out_dir <- paste0(gene, "/figure/", gene, "/", intron, "@", SNP)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- paste0(out_dir, "/", gene, "@", intron, "@", SNP, ".pdf")
ggsave(out_pdf, wrap_plots(plots, ncol = 4), width = 18, height = 12)


