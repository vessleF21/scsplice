library(data.table)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
SNP <- args[1]
gene <- args[2]
gene_id <- args[3]

celltypes <- c("BNav", "BMem", "CD4EM", "CD4Naive", "CD4Treg", "CD8Naive", "CD8GZMH", "CD8GZMK", "MAIT", "NKDim", "NKBright", "MonocM", "NoClaM")

color_mapping <- c(
  "BNav" = "#01579B", "BMem" = "#A5EDFF", "CD4EM" = "#1B5E20", "CD4Naive" = "#8BC34A",
  "CD4Treg" = "#D4E157", "CD8GZMH" = "#FF5722", "CD8GZMK" = "#FFCDD2", "CD8Naive" = "#F9A825",
  "MAIT" = "#FFEE58", "MonocM" = "#6A1B9A", "NoClaM" = "#E1BEE7", "NKBright" = "#82581F", "NKDim" = "#D7CCC8"
)

genocol <- read.table(paste0(gene, "/snp/snp_head.txt"), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
genotype <- read.table(paste0(gene, "/snp/", SNP, ".txt"), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
colnames(genotype) <- genocol
REF <- genotype$REF[1]
ALT <- genotype$ALT[1]

genotype[10:ncol(genotype)] <- lapply(genotype[10:ncol(genotype)], function(col) sapply(col, function(gt) {
  if (gt %in% c("0/1", "1/0")) paste0(REF, ALT)
  else if (gt == "0/0") paste0(REF, REF)
  else if (gt == "1/1") paste0(ALT, ALT)
  else NA
}))

plots <- list()

for (celltype in celltypes) {
  color <- color_mapping[[celltype]]
  file_path <- paste0(gene, "/phe/data/", celltype, "_", gene_id, ".txt")
  
  if (!file.exists(file_path)) {
    p <- ggplot() +
      labs(title = paste("No Data for", celltype), x = "Genotype", y = "Phenotype") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    plots[[length(plots) + 1]] <- p
    next
  }
  
  col_phe <- as.character(unlist(read.table(paste0(gene, "/phe/data/", celltype, "_head.txt"), sep = "\t", check.names = FALSE)))
  phenotype <- fread(file_path, sep = "\t")
  colnames(phenotype) <- col_phe
  
  genotype <- as.data.frame(genotype)
  phenotype <- as.data.frame(phenotype)
  geno_phe <- intersect(as.character(colnames(genotype)), as.character(colnames(phenotype)))
  merged_data <- rbind(
    phenotype[, geno_phe, drop = FALSE],
    genotype[, geno_phe, drop = FALSE]
  )
  
  transposed_data <- as.data.frame(t(merged_data))
  rownames(transposed_data) <- NULL
  colnames(transposed_data) <- c("phe", "geno")
  transposed_data$geno <- factor(transposed_data$geno, levels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)))
  transposed_data$phe <- as.numeric(transposed_data$phe)
  
  p <- ggplot(transposed_data, aes(x = geno, y = phe)) +
    geom_violin(trim = FALSE, alpha = 0, scale = "width", size = 0.5, color = color) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = color) +
    labs(title = paste(celltype), x = "Genotype", y = "Phenotype") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      axis.line.x = element_line(color = "black", size = 0.5),
      axis.line.y = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      panel.grid = element_blank()
    )
  
  plots[[length(plots) + 1]] <- p
}

output_dir <- paste0(gene, "/figure/", gene, "/", celltype, "@", SNP)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- paste0(output_dir, "/", gene, "@", celltype, "@", SNP, ".pdf")
combined_plot <- wrap_plots(plots, ncol = 4)
ggsave(output_file, combined_plot, width = 18, height = 12)
