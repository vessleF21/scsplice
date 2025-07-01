color_mapping <- c(
  "BNav" = "#01579B",
  "BMem" = "#A5EDFF",
  "CD4EM" = "#1B5E20",
  "CD4Naive" = "#8BC34A",
  "CD4Treg" = "#D4E157",
  "CD8GZMH" = "#FF5722",
  "CD8GZMK" = "#FFCDD2",
  "CD8Naive" = "#F9A825",
  "MAIT" = "#FFEE58",
  "MonocM" = "#6A1B9A",
  "NoClaM" = "#E1BEE7",
  "NKBright" = "#82581F",
  "NKDim" = "#D7CCC8"
)

args <- commandArgs(trailingOnly = TRUE)
SNP <- args[1]
celltype <- args[2]
color <- color_mapping[[celltype]]
intron <- as.character(args[3])
antype <- args[4]
gene <- args[5]
chr = as.character(sub("chr(\\d+):.*", "\\1", SNP))

library(data.table)
library(ggplot2)
library(patchwork)

genocol <- read.table(paste0("geno.txt"), sep = "\t", check.names = FALSE)
genotype <- read.table(paste0("SNP_info/", SNP, ".txt"), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
colnames(genotype) <- genocol
REF <- genotype$REF[1]
ALT <- genotype$ALT[1]
genotype[10:ncol(genotype)] <- lapply(genotype[10:ncol(genotype)], function(col) sapply(col, function(gt) if (gt %in% c("0/1", "1/0")) paste0(REF, ALT) else if (gt == "0/0") paste0(REF, REF) else if (gt == "1/1") paste0(ALT, ALT) else NA))

col_phe <- read.table(paste("phename/", celltype, "_phename.txt", sep=""), sep=" ", check.names=FALSE)
newphe <- as.data.frame(matrix(NA,0,ncol(col_phe)))
phenotype <- fread(paste("S20_", celltype, "/S20_", celltype, chr, "_phenotype.txt.gz", sep=""), sep="\t")
phenotype <- as.data.frame(phenotype)
phenotype <- phenotype[phenotype$ID == intron, ]

PC <- read.table(paste("S20_", celltype, "/", celltype, "_ok_PC.txt", sep=""), sep="\t", header=TRUE, row.names=1, check.names=FALSE)
PC <- PC[c(14),]
PC_low45 <- colnames(PC)[as.numeric(PC[1, ]) < 45]
PC_M4565 <- colnames(PC)[as.numeric(PC[1, ]) >= 45 & as.numeric(PC[1, ]) <= 65]
PC_more65 <- colnames(PC)[as.numeric(PC[1, ]) > 65]

H_genocol <- intersect(as.character(colnames(genotype)), as.character(PC_more65))
H_geno <- genotype[, colnames(genotype) %in% H_genocol]

M_genocol <- intersect(as.character(colnames(genotype)), as.character(PC_M4565))
M_geno <- genotype[, colnames(genotype) %in% M_genocol]

Y_genocol <- intersect(as.character(colnames(genotype)), as.character(PC_low45))
Y_geno <- genotype[, colnames(genotype) %in% Y_genocol]

H_phe <- phenotype[, colnames(phenotype) %in% H_genocol]
M_phe <- phenotype[, colnames(phenotype) %in% M_genocol]
Y_phe <- phenotype[, colnames(phenotype) %in% Y_genocol]

H_common_columns <- intersect(colnames(H_geno), colnames(H_phe))
H_merged_data <- rbind(H_phe[, H_common_columns, drop = FALSE], H_geno[, H_common_columns, drop = FALSE])
H_transposed_data <- as.data.frame(t(H_merged_data))
rownames(H_transposed_data) <- NULL
colnames(H_transposed_data) <- c("phe","geno")
H_transposed_data$geno <- factor(H_transposed_data$geno, levels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)), labels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)))
H_transposed_data$phe <- as.numeric(H_transposed_data$phe)

M_common_columns <- intersect(colnames(M_geno), colnames(M_phe))
M_merged_data <- rbind(M_phe[, M_common_columns, drop = FALSE], M_geno[, M_common_columns, drop = FALSE])
M_transposed_data <- as.data.frame(t(M_merged_data))
rownames(M_transposed_data) <- NULL
colnames(M_transposed_data) <- c("phe","geno")
M_transposed_data$geno <- factor(M_transposed_data$geno, levels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)), labels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)))
M_transposed_data$phe <- as.numeric(M_transposed_data$phe)

Y_common_columns <- intersect(colnames(Y_geno), colnames(Y_phe))
Y_merged_data <- rbind(Y_phe[, Y_common_columns, drop = FALSE], Y_geno[, Y_common_columns, drop = FALSE])
Y_transposed_data <- as.data.frame(t(Y_merged_data))
rownames(Y_transposed_data) <- NULL
colnames(Y_transposed_data) <- c("phe","geno")
Y_transposed_data$geno <- factor(Y_transposed_data$geno, levels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)), labels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)))
Y_transposed_data$phe <- as.numeric(Y_transposed_data$phe)

output_dir <- paste0(intron)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
setwd(output_dir)


H_p <- ggplot(H_transposed_data, aes(x = geno, y = phe)) +
  geom_violin(trim = FALSE, alpha = 0, scale = "width", size = 0.5, color = color) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = color) +
  labs(title = paste("Violin Plot of", celltype), x = "Genotype", y = "Phenotype") +
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

M_p <- ggplot(M_transposed_data, aes(x = geno, y = phe)) +
  geom_violin(trim = FALSE, alpha = 0, scale = "width", size = 0.5, color = color) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = color) +
  labs(title = paste("Violin Plot of", celltype), x = "Genotype", y = "Phenotype") +
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

Y_p <- ggplot(Y_transposed_data, aes(x = geno, y = phe)) +
  geom_violin(trim = FALSE, alpha = 0, scale = "width", size = 0.5, color = color) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = color) +
  labs(title = paste("Violin Plot of", celltype), x = "Genotype", y = "Phenotype") +
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

combined_plot <- H_p | M_p | Y_p
ggsave("combined_plot.pdf", plot = combined_plot, width = 18, height = 6, dpi = 300)
