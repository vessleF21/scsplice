

library(ggplot2)
library(data.table)
library(dplyr)
library(reshape2)

GENE        <- "CTLA4"
CLUSTER     <- "chr2:clu_12193_+"
INPUT_FILE  <- "oneK1K_perind.counts_merged_with_cluster.txt.marker_ok"
CELL_REF    <- "cell_type_use.txt"
OUTPUT_CSV  <- "dat_melt_group_tmp_dcast_tmp.csv"
OUTPUT_PDF  <- "Fig.2b_CTLA4.pdf"

IDS_TO_MATCH <- c(
  "chr2:203870933:203872708",
  "chr2:203871487:203872708"
)

IID_LABEL <- c(
  "CTLA4__chr2:203870933:203872708:clu_12193_+" = "CTLA4: exon3-exon5 (sCTLA4)",
  "CTLA4__chr2:203871487:203872708:clu_12193_+" = "CTLA4: exon4-exon5 (CTLA4)"
)

ORD_LEVEL_CELL <- c("CD4Naive", "CD4Treg", "CD4EM")

ORD_LEVEL_IID <- c(
  "CTLA4: exon3-exon5 (sCTLA4)",
  "CTLA4: exon4-exon5 (CTLA4)"
)

celltype_pal <- c(
  "BMem"     = "#FF5733",
  "BNav"     = "#33FF57",
  "CD4EM"    = "#3357FF",
  "CD4Naive" = "#FF33FF",
  "CD4Treg"  = "#FF33A1",
  "CD8GZMH"  = "#33FFA1",
  "CD8GZMK"  = "#A1FF33",
  "CD8Naive" = "#FFA133",
  "MAIT"     = "#5733FF",
  "MonocM"   = "#A133FF",
  "NKBright" = "#33A1FF",
  "NKDim"    = "#FF5733",
  "NoClaM"   = "#33FF57"
)

celltype_colors <- c(
  "CD4EM"    = "#1B5E20",
  "CD4Naive" = "#8BC34A",
  "CD4Treg"  = "#D4E157"
)

label_dict <- list(
  "BMem"     = "BMem",
  "BNav"     = "BNav",
  "CD4EM"    = "CD4EM",
  "CD4Naive" = "CD4Naive",
  "CD4Treg"  = "CD4Treg",
  "CD8GZMH"  = "CD8GZMH",
  "CD8GZMK"  = "CD8GZMK",
  "CD8Naive" = "CD8Naive",
  "MAIT"     = "MAIT",
  "MonocM"   = "MonocM",
  "NKBright" = "NKBright",
  "NKDim"    = "NKDim",
  "NoClaM"   = "NoClaM"
)

celltypes <- read.table(CELL_REF)$V1

dat <- read.delim(INPUT_FILE, sep = " ", check.names = FALSE)

dat.melt <- melt(dat,
                 id.vars    = c(1, 2, ncol(dat)),
                 value.name = "normalized_splicing")

colnames(dat.melt)[2] <- "ID"

dat.melt$celltype <- gsub(
  "\\.bam$", "",
  sapply(strsplit(as.character(dat.melt$variable), "@"),
         function(x) x[length(x)])
)

dat.melt$normalized_splicing <- sapply(
  dat.melt[, 5],
  function(x) eval(parse(text = x))
)

for (i in seq_along(label_dict)) {
  pat <- paste0("^", names(label_dict)[i], "$")
  dat.melt$celltype   <- gsub(pat, label_dict[[i]], dat.melt$celltype)
  names(celltype_pal) <- gsub(pat, label_dict[[i]], names(celltype_pal))
}

dat.melt <- na.omit(dat.melt)

dat.melt.group <- dat.melt |>
  dplyr::group_by(ID, celltype, genes) |>
  dplyr::summarise(
    average_normalized_splicing = mean(normalized_splicing),
    .groups = "drop"
  ) |>
  dplyr::mutate(iid = paste(genes, ID, sep = "__"))

CTLA4 <- dat.melt.group |>
  dplyr::filter(genes == GENE) |>
  dplyr::filter(grepl(paste(IDS_TO_MATCH, collapse = "|"), ID))


dat.tmp <- CTLA4 |>
  dplyr::filter(
    celltype %in% grep("CD4|CD8|MAIT",
                       unique(dat.melt.group$celltype), value = TRUE)
  ) |>
  dplyr::mutate(
    iid_label = IID_LABEL[iid],
    celltype  = factor(celltype, levels = ORD_LEVEL_CELL),
    iid_label = factor(iid_label, levels = ORD_LEVEL_IID)
  )

dat.wide <- dcast(dat.tmp,
                  celltype ~ iid_label,
                  value.var = "average_normalized_splicing")
dat.wide[is.na(dat.wide)] <- 0

dat.plot <- melt(dat.wide,
                 id.vars       = 1,
                 variable.name = "iid_label",
                 value.name    = "average_normalized_splicing") |>
  dplyr::filter(iid_label != "CTLA4: exon3-exon4 (CTLA4)") |>
  dplyr::mutate(iid_label = factor(iid_label, levels = ORD_LEVEL_IID))

write.csv(dat.plot, file = OUTPUT_CSV, row.names = FALSE)

f1 <- ggplot(dat.plot,
             aes(x     = iid_label,
                 y     = average_normalized_splicing,
                 color = celltype)) +
  geom_point(size = 4) +
  scale_color_manual(values = celltype_colors) +
  ylim(0, 0.85) +
  labs(
    x     = "IID Label",
    y     = "Average Normalized Splicing",
    color = "Cell Type"
  ) +
  theme_bw() +
  theme(
    legend.position  = "top",
    legend.title     = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y      = element_text(size = 10),
    axis.title.x     = element_text(size = 12, face = "bold"),
    axis.title.y     = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank()
  )

print(f1)
ggsave(OUTPUT_PDF, f1, width = 6, height = 6)