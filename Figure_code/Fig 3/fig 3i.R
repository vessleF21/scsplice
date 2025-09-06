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

genocol <- read.table(paste0(,antype,"/geno.txt",sep=""), sep = "\t", check.names = FALSE)
genotype <- read.table(paste0(,antype,"/SNP_info/", SNP, ".txt",sep=""), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")

colnames(genotype)<-genocol
REF <- genotype$REF[1]
ALT <- genotype$ALT[1]
genotype[10:ncol(genotype)] <- lapply(genotype[10:ncol(genotype)], function(col) sapply(col, function(gt) if (gt %in% c("0/1", "1/0")) paste0(REF, ALT) else if (gt == "0/0") paste0(REF, REF) else if (gt == "1/1") paste0(ALT, ALT) else NA))

col_phe<-read.table(paste(antype,"/phename/",celltype,"_phename.txt",sep=""),sep=" ",check.names=F)
newphe<-as.data.frame(matrix(NA,0,ncol(col_phe)))
phenotype<-fread(paste(celltype,"/S20_",celltype,chr,"_phenotype.txt.gz",sep=""),sep="\t")
phenotype<-as.data.frame(phenotype)

phenotype <- phenotype[phenotype$ID == intron, ]

PC<-read.table(paste(antype,"/PC/",celltype,"_ok_PC.txt",sep=""),sep="\t",header=TRUE,row.names=1,check.names=F)

PC<-PC[c(13),]
fe_PC<-colnames(PC)[PC[1, ] == 1] 
male_PC<-colnames(PC)[PC[1, ] == 0] 

fe_genocol <- intersect(as.character(colnames(genotype)), as.character(fe_PC))
fe_geno <- genotype[, colnames(genotype) %in% fe_genocol]
male_genocol <- intersect(as.character(colnames(genotype)), as.character(male_PC))
male_geno <- genotype[, colnames(genotype) %in% male_genocol]
dim(male_geno)
dim(fe_geno)

fe_phe <- phenotype[, colnames(phenotype) %in% fe_genocol]
male_phe <- phenotype[, colnames(phenotype) %in% male_genocol]
dim(fe_phe)
dim(male_phe)

male_common_columns <- intersect(colnames(male_geno), colnames(male_phe))
male_merged_data <- rbind(male_phe[, male_common_columns, drop = FALSE],male_geno[, male_common_columns, drop = FALSE])
male_transposed_data <- as.data.frame(t(male_merged_data))
rownames(male_transposed_data) <- NULL
colnames(male_transposed_data) <- c("phe","geno")

male_transposed_data$geno <- factor(male_transposed_data$geno, levels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)), labels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)))
male_transposed_data$phe <- as.numeric(male_transposed_data$phe)

fe_common_columns <- intersect(colnames(fe_geno), colnames(fe_phe))
fe_merged_data <- rbind(fe_phe[, fe_common_columns, drop = FALSE],fe_geno[, fe_common_columns, drop = FALSE])
fe_transposed_data <- as.data.frame(t(fe_merged_data))
rownames(fe_transposed_data) <- NULL
colnames(fe_transposed_data) <- c("phe","geno")

fe_transposed_data$geno <- factor(fe_transposed_data$geno, levels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)), labels = c(paste0(REF, REF), paste0(REF, ALT), paste0(ALT, ALT)))
fe_transposed_data$phe <- as.numeric(fe_transposed_data$phe)

output_dir <- paste0(intron)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
setwd(output_dir)

library(dplyr)
male_p <- ggplot(male_transposed_data, aes(x = geno, y = phe)) +
  geom_violin(trim = FALSE, alpha = 0, scale = "width", size = 0.5, color = color) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = color) +
  labs(title = paste("Violin Plot of", celltype),
       x = "Genotype",
       y = "Phenotype") +
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
