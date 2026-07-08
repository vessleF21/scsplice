

library(data.table)
library(ggplot2)
library(patchwork)

args   <- commandArgs(trailingOnly = TRUE)
SNP    <- args[1]
intron <- args[2]
gene   <- args[3]

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

path_snp     <- function(gene, SNP)         paste0(gene, "/snp/", SNP, ".txt")
path_phename <- function(ct)                paste0("phename/", ct, "_phename.txt")
path_phe     <- function(gene, intron, SNP, ct) paste0(gene, "/phe/cellty/", gene, "/", intron, "@", SNP, "/", ct, "@", intron, "@", SNP, ".txt")
path_out_dir <- function(gene, intron, SNP)     paste0(gene, "/figure/", gene, "/", intron, "@", SNP)

genocol  <- read.table("geno.txt", sep = "\t", check.names = FALSE)
genotype <- read.table(path_snp(gene, SNP),
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

plots               <- list()
sample_summary_list <- list()

for (ct in CELL_ORDER) {
  color     <- COLORS_CELL[[ct]]
  file_path <- path_phe(gene, intron, SNP, ct)

  if (!file.exists(file_path)) {
    cat(sprintf("[%s] file not found\n", ct))
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
    cat(sprintf("[%s] no common samples\n", ct))
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
  cat(sprintf("[%s]  total: %d\n", ct, sum(sizes)))
  print(sizes)

  for (geno_lvl in names(sizes)) {
    sample_summary_list[[length(sample_summary_list) + 1]] <-
      data.frame(CellType = ct, Genotype = geno_lvl,
                 n = as.numeric(sizes[geno_lvl]), Status = "OK")
  }

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

sample_summary <- do.call(rbind, sample_summary_list)

celltype_texts <- sapply(CELL_ORDER, function(ct) {
  ok <- sample_summary[sample_summary$CellType == ct & sample_summary$Status == "OK", ]
  if (nrow(ok) > 0) {
    paste0(ct, " - ", paste(ok$Genotype, " (n=", ok$n, ")", sep = "", collapse = ", "))
  } else {
    st <- unique(sample_summary$Status[sample_summary$CellType == ct])
    paste0(ct, " - ", ifelse(length(st) > 0 && st[1] == "File not found", "File not found", "No data"))
  }
})

total_n <- sum(sample_summary$n[sample_summary$Status == "OK"], na.rm = TRUE)
combined_text <- paste0("Sample sizes by cell type and genotype: ",
                        paste(celltype_texts, collapse = "; "),
                        ". Total n = ", total_n)

cat(celltype_texts, sep = "\n")
cat(sprintf("\nTotal n: %d\n", total_n))
cat(combined_text, "\n")

out_dir <- path_out_dir(gene, intron, SNP)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_summary$SNP    <- SNP
sample_summary$Gene   <- gene
sample_summary$Intron <- intron
sample_summary <- sample_summary[, c("SNP", "Gene", "Intron", "CellType", "Genotype", "n", "Status")]

write.table(sample_summary, file.path(out_dir, "sample_sizes.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

writeLines(c(
  paste("SNP:", SNP), paste("Gene:", gene),
  paste("Intron:", intron), paste("REF:", REF, "ALT:", ALT),
  "", celltype_texts, "",
  paste("Total n =", total_n), "", combined_text
), file.path(out_dir, "sample_sizes_formatted.txt"))

out_pdf <- file.path(out_dir, paste0(gene, "@", intron, "@", SNP, ".pdf"))
ggsave(out_pdf, wrap_plots(plots, ncol = 4), width = 18, height = 12)

message("saved: ", out_pdf)


###run
Rscript plot_script.R <SNP> <intron> <gene>

Rscript plot_script.R chr2:230251261:C:T chr2:230245080:230247916_+ SP140


