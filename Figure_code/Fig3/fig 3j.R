


library(dplyr)
library(data.table)

# args <- commandArgs(trailingOnly = TRUE)
# gwas_name <- args[1]; celltype <- args[2]; chr <- args[3]
# intron <- args[4]; fenc1 <- args[5]; fenc2 <- args[6]

gwas_name <- "ChenMH_CELL_2020_BAS"
celltype  <- "NoClaM"
chr       <- "11"
intron    <- "chr11:319985:320565:clu_84_-"
fenc1     <- "male"
fenc2     <- "female"

snp_1mb_file   <- paste0("<path>/", celltype, "_", intron, ".txt")
snp_1mb_file_M <- paste0("<path>/", fenc1, "_", celltype, "_", intron, ".txt")
snp_1mb_file_F <- paste0("<path>/", fenc2, "_", celltype, "_", intron, ".txt")
GWAS_data      <- paste0("<path>/", gwas_name, "/", chr, ".txt.gz")
maf_file       <- paste0("<path>/maf_chr", chr, ".txt")
id_tr          <- paste0("<path>/id_tr_rs_chr", chr, ".txt")

print("Start reading input data...")

snp_1mb <- read.table(snp_1mb_file, header = FALSE, sep = " ", stringsAsFactors = FALSE)
colnames(snp_1mb) <- c("gene", "intron", "variant_id", "pval_nominal", "beta")
snp_1mb$variant_id <- gsub("^chr", "", snp_1mb$variant_id)

GWAS_pre <- fread(GWAS_data, sep = "\t", header = TRUE) |> as.data.frame()
GWAS_pre <- GWAS_pre[, 2:5]
colnames(GWAS_pre) <- c("variant_id", "beta", "se", "pval_nominal")

snp_1mb_M <- fread(snp_1mb_file_M, sep = " ", header = FALSE) |> as.data.frame()
snp_1mb_M <- snp_1mb_M[, c(3, 4)]
colnames(snp_1mb_M) <- c("variant_id", "pval_nominal_M")
snp_1mb_M$variant_id <- gsub("^chr", "", snp_1mb_M$variant_id)

snp_1mb_F <- fread(snp_1mb_file_F, sep = " ", header = FALSE) |> as.data.frame()
snp_1mb_F <- snp_1mb_F[, c(3, 4)]
colnames(snp_1mb_F) <- c("variant_id", "pval_nominal_F")
snp_1mb_F$variant_id <- gsub("^chr", "", snp_1mb_F$variant_id)

MAF <- fread(maf_file, header = TRUE) |> as.data.frame()
colnames(MAF) <- c("variant_id", "maf")

id_tr_rs <- read.table(id_tr, header = FALSE, sep = " ", stringsAsFactors = FALSE)
colnames(id_tr_rs) <- c("variant_id_new", "rs")

merged <- snp_1mb |>
  merge(MAF,      by = "variant_id") |>
  merge(GWAS_pre, by = "variant_id", suffixes = c("_sqtl", "_gwas")) |>
  dplyr::select(variant_id, pval_nominal_sqtl, pval_nominal_gwas) |>
  merge(snp_1mb_M, by = "variant_id") |>
  merge(snp_1mb_F, by = "variant_id") |>
  dplyr::mutate(
    variant_id_new = sapply(strsplit(variant_id, ":", fixed = TRUE),
                            \(x) paste0(x[1], ":", x[2]))
  ) |>
  merge(id_tr_rs, by = "variant_id_new")

out_dir <- paste(gwas_name, celltype, intron, sep = "@")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

.save <- function(dat, col, fname) {
  out <- dat |> dplyr::select(rs, P = {{ col }})
  write.table(out, file.path(out_dir, fname),
              row.names = FALSE, sep = "\t", quote = FALSE)
}

.save(merged, pval_nominal_gwas, "GWAS_data.txt")
.save(merged, pval_nominal_sqtl, "all_sqtl.txt")
.save(merged, pval_nominal_M,    "M_sqtl.txt")
.save(merged, pval_nominal_F,    "F_sqtl.txt")




This figure was generated using the LocusZoom web tool (http://locuszoom.org/).


