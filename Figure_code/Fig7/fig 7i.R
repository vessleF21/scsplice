

rm(list = ls())
gc()

library(dplyr)
library(ggplot2)
library(data.table)
library(locuscomparer)

BASE_DIR    <- "IVNS1ABP_RPS24"
RSID_FILE   <- file.path(BASE_DIR, "rsid_ref.txt")
MERGED_FILE <- file.path(BASE_DIR, "merged_IVNS1ABP_RPS24_184296388_186317273.txt")
LD_FILE     <- file.path(BASE_DIR, "results/ld.tsv")

LABEL1 <- "IVNS1ABP cis-eQTL"
LABEL2 <- "RPS24 trans-sQTL"
CHR    <- 1
SNP_ID <- "rs10798014"

total  <- fread(RSID_FILE)
merged <- fread(MERGED_FILE, sep = "\t")
ld     <- fread(LD_FILE, sep = "\t")

merged <- merge(total, merged, by.x = "Variation ID", by.y = "VID", all = FALSE)

merged <- merged %>%
  rename(rsid = dbSNP) %>%
  select(rsid, pos, epval, spval, logpe, logps)
colnames(merged) <- c("rsid", "pos", "pval1", "pval2", "logp1", "logp2")

merged$chr <- as.character(CHR)
merged     <- merged[!duplicated(merged$pos), ]

snp   <- merged[which.min(merged$pval1 * merged$pval2), "rsid"][[1]]
snp   <- get_lead_snp(merged, snp)
color <- assign_color(merged$rsid, snp, ld)
shape <- setNames(ifelse(merged$rsid == snp, 23, 21), merged$rsid)
size  <- setNames(ifelse(merged$rsid == snp,  4,  3), merged$rsid)
merged$label <- ifelse(merged$rsid == snp, merged$rsid, "")

legend          <- TRUE
legend_position <- c("bottomright", "topright", "topleft")

p_scatter <- make_scatterplot(merged, LABEL1, LABEL2,
                               color, shape, size, legend, legend_position)
ggsave("single.pdf", p_scatter, width = 4, height = 4)

p2 <- make_combined_plot(merged, LABEL1, LABEL2,
                          ld, chr = CHR, snp = SNP_ID,
                          combine = FALSE, legend)
ggsave("locuscompare.pdf", plot = p2$locuscompare, width = 6, height = 6)
ggsave("locuszoom1.pdf",   plot = p2$locuszoom1,   width = 6, height = 6)
ggsave("locuszoom2.pdf",   plot = p2$locuszoom2,   width = 6, height = 6)

