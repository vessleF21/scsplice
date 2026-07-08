

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

path_geno    <- function()                      "geno.txt"
path_snp     <- function(gene, SNP)             paste0("/ggvilllin/", gene, "/snp/", SNP, ".txt")
path_phename <- function(ct)                    paste0("/phename/", ct, "_phename.txt")
path_phe     <- function(gene, intron, SNP, ct) paste0("/ggvilllin/", gene, "/phe/cellty/", gene, "/", intron, "@", SNP, "/", ct, "@", intron, "@", SNP, ".txt")
path_out_dir <- function(gene, intron, SNP)     paste0("/ggvilllin/", gene, "/figure/", gene, "/", intron, "@", SNP)

args   <- commandArgs(trailingOnly = TRUE)
SNP    <- args[1]
intron <- args[2]
gene   <- args[3]

genocol  <- read.table(path_geno(), sep = "\t", check.names = FALSE)
genotype <- read.table(path_snp(gene, SNP), sep = "\t", check.names = FALSE,
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

plots               <- list()
sample_summary_list <- list()
desc_stats_list     <- list()

for (ct in CELL_ORDER) {
  color     <- COLORS_CELL[[ct]]
  file_path <- path_phe(gene, intron, SNP, ct)

  if (!file.exists(file_path)) {
    plots[[length(plots) + 1]] <- ggplot() +
      labs(title = paste("No Data:", ct), x = "Genotype", y = "Phenotype") +
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
    sample_summary_list[[length(sample_summary_list) + 1]] <-
      data.frame(CellType = ct, Genotype = NA, n = 0, Status = "File not found")
    next
  }

  col_phe   <- as.character(unlist(read.table(path_phename(ct), sep = " ", check.names = FALSE)))
  phenotype <- as.data.frame(fread(file_path, sep = " "))
  colnames(phenotype) <- col_phe
  genotype  <- as.data.frame(genotype)

  common_samples <- intersect(colnames(genotype), colnames(phenotype))[-1]

  if (length(common_samples) == 0) {
    plots[[length(plots) + 1]] <- ggplot() +
      labs(title = paste("No Common Samples:", ct), x = "Genotype", y = "Phenotype") +
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
    sample_summary_list[[length(sample_summary_list) + 1]] <-
      data.frame(CellType = ct, Genotype = NA, n = 0, Status = "No common samples")
    next
  }

  merged <- as.data.frame(t(rbind(
    phenotype[, common_samples, drop = FALSE],
    genotype[,  common_samples, drop = FALSE]
  )))
  rownames(merged) <- NULL
  colnames(merged) <- c("phe", "geno")
  merged$geno <- factor(merged$geno,
                        levels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)))
  merged$phe  <- as.numeric(merged$phe)

  sizes <- table(merged$geno, useNA = "ifany")
  for (geno_lvl in names(sizes)) {
    sample_summary_list[[length(sample_summary_list) + 1]] <-
      data.frame(CellType = ct, Genotype = geno_lvl,
                 n = as.numeric(sizes[geno_lvl]), Status = "OK")
  }

  for (lvl in levels(merged$geno)) {
    vals <- merged$phe[merged$geno == lvl & !is.na(merged$geno)]
    if (length(vals) == 0) next
    desc_stats_list[[length(desc_stats_list) + 1]] <- data.frame(
      CellType = ct, Genotype = lvl, n = length(vals),
      Mean     = mean(vals, na.rm = TRUE),
      Median   = median(vals, na.rm = TRUE),
      SD       = sd(vals, na.rm = TRUE),
      Min      = min(vals, na.rm = TRUE),
      Max      = max(vals, na.rm = TRUE),
      Q25      = quantile(vals, 0.25, na.rm = TRUE),
      Q75      = quantile(vals, 0.75, na.rm = TRUE)
    )
  }

  plots[[length(plots) + 1]] <- ggplot(merged, aes(x = geno, y = phe)) +
    geom_violin(trim = FALSE, alpha = 0, scale = "width",
                linewidth = 0.5, color = color) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = color) +
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

sample_summary <- do.call(rbind, sample_summary_list)
desc_stats     <- do.call(rbind, desc_stats_list)

add_meta <- function(df) { df$SNP <- SNP; df$Gene <- gene; df$Intron <- intron; df }
sample_summary <- add_meta(sample_summary)[, c("SNP","Gene","Intron","CellType","Genotype","n","Status")]
desc_stats     <- add_meta(desc_stats)[,    c("SNP","Gene","Intron","CellType","Genotype","n","Mean","Median","SD","Min","Max","Q25","Q75")]

total_n <- sum(sample_summary$n[sample_summary$Status == "OK"], na.rm = TRUE)

make_texts <- function(ct, with_stats = FALSE) {
  ok    <- sample_summary[sample_summary$CellType == ct & sample_summary$Status == "OK", ]
  stats <- desc_stats[desc_stats$CellType == ct, ]
  if (nrow(ok) == 0) {
    st <- unique(sample_summary$Status[sample_summary$CellType == ct])
    return(paste0(ct, " - ", ifelse(length(st) > 0 && st[1] == "File not found", "File not found", "No data")))
  }
  if (with_stats && nrow(stats) > 0) {
    paste0(ct, " - ", paste(sprintf("%s (n=%d, mean=%.3f, median=%.3f)",
                                    stats$Genotype, stats$n, stats$Mean, stats$Median), collapse = ", "))
  } else {
    paste0(ct, " - ", paste(ok$Genotype, " (n=", ok$n, ")", sep = "", collapse = ", "))
  }
}

texts_simple   <- sapply(CELL_ORDER, make_texts, with_stats = FALSE)
texts_detailed <- sapply(CELL_ORDER, make_texts, with_stats = TRUE)

out_dir <- path_out_dir(gene, intron, SNP)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.table(sample_summary, file.path(out_dir, "sample_sizes.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(desc_stats, file.path(out_dir, "descriptive_stats.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

writeLines(c(
  paste("SNP:", SNP), paste("Gene:", gene),
  paste("Intron:", intron), paste("REF:", REF, "ALT:", ALT),
  "", "=== Sample Sizes ===", texts_simple, "", paste("Total n =", total_n),
  "", "=== Descriptive Statistics ===", texts_detailed,
  "", paste0("Sample sizes by cell type and genotype: ",
             paste(texts_simple, collapse = "; "), ". Total n = ", total_n),
  "", paste0("Descriptive statistics by cell type and genotype: ",
             paste(texts_detailed, collapse = "; "), ". Total n = ", total_n)
), file.path(out_dir, "summary.txt"))

out_pdf <- file.path(out_dir, paste0(gene, "@", intron, "@", SNP, ".pdf"))
ggsave(out_pdf, wrap_plots(plots, ncol = 4), width = 18, height = 12)

