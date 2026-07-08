

grep SYNRG $GTF > SYNRG.gencode.txt

python3 - << 'PYEOF'
import re
from collections import defaultdict
import pandas as pd

INPUT  = "SYNRG.gencode.txt"
OUTPUT = "SYNRG_model.tsv"

records = defaultdict(list)
with open(INPUT) as f:
    for line in f:
        cols      = line.split("\t")
        feat_type = cols[2]
        if feat_type != "exon":
            continue
        annotation = cols[8]
        transcript = re.search(r"ENST[0-9]+", annotation).group(0)
        records["chrom"].append(cols[0])
        records["transcript"].append(transcript)
        records["start"].append(cols[3])
        records["end"].append(cols[4])
        records["strand"].append(cols[6])
        records["orientation"].append(1)

pd.DataFrame(records).to_csv(OUTPUT, sep="\t", index=False)
PYEOF

awk '$2 == "transcript" || $2 == "ENST00000612223" || $2 == "ENST00000621136"' \
    SYNRG_model.tsv > SYNRG_newmodel.tsv

cat > SYNRG_features.tsv << 'TSVEOF'
transcript  name  type  position  forward
ENST00000612223 tss tss 37596344  FALSE
ENST00000621136 tss tss 37596344  FALSE
TSVEOF

Rscript - << 'REOF'
library(ggplot2)
library(gggenes)
library(data.table)
library(dplyr)

SNP_POS  <- 37596344L
OUT_FILE <- "SYNRG_model_z.pdf"

tchp     <- fread("SYNRG_newmodel.tsv")
features <- fread("SYNRG_features.tsv")

tchp <- tchp %>%
  mutate(is_target_exon = (start == SNP_POS))

p <- ggplot(tchp, aes(xmin = start, xmax = end, y = transcript,
                      fill = is_target_exon, color = is_target_exon)) +
  geom_gene_arrow(
    arrowhead_width  = grid::unit(0, "mm"),
    arrowhead_height = grid::unit(3, "mm")
  ) +
  scale_fill_manual(values  = c("FALSE" = "black", "TRUE" = "red")) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_feature(
    data = features,
    aes(x = position, y = transcript, forward = forward)
  ) +
  theme_genes() +
  theme(legend.position = "none") +
  labs(x = "chr17", y = NULL)

ggsave(OUT_FILE, p, width = 6, height = 3)
REOF



