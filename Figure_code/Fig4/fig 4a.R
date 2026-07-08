

BASE="test"
OUT="${BASE}/nc2/1/strength"
GENE_LIST="${OUT}/gene_list.txt"
NOM_DIR="${BASE}/splic_split/junc"
CELL_ORDER="BNav BMem CD4Naive CD4EM CD4Treg CD8Naive CD8GZMH CD8GZMK MAIT NKDim NKBright MonocM NoClaM"

mkdir -p "${OUT}" && cd "${OUT}" || exit 1

# Step 1: generate per-gene beta extraction script
cat > script.sh << EOF
#!/bin/bash
cluster="\$1"; snp="\$2"; gene="\$3"; chrnum="\$4"

grep -H "\${cluster}" ${NOM_DIR}/S20_*/nominal7/*_"\${chrnum}"_* \\
| grep "\${snp}" \\
| grep "\${gene}" \\
| awk -v order="${CELL_ORDER}" '
    BEGIN {
        n = split(order, ord, " ")
        for (i = 1; i <= n; i++) { want[ord[i]] = 1; val[ord[i]] = 0 }
    }
    {
        fname = \$1; sub(/^.*\//, "", fname)
        split(fname, a, "_"); sample = a[1]
        if (sample in want) val[sample] = \$(NF-1)
    }
    END {
        for (i = 1; i <= n; i++) print ord[i], val[ord[i]]
    }' > "\${gene}.txt"
EOF
chmod +x script.sh

# Step 2: extract beta values per gene
# gene_list.txt format: <cluster> <snp> <gene> <chr>
while read -r cluster snp gene chrnum; do
    [ -z "${cluster}" ] && continue
    bash script.sh "${cluster}" "${snp}" "${gene}" "${chrnum}"
done < "${GENE_LIST}"

# Step 3: merge beta values into matrix (rows=cell types, cols=genes)
mapfile -t genes < <(awk '{print $3}' "${GENE_LIST}")

col_files=()
for g in "${genes[@]}"; do
    [ -f "${g}.txt" ] || { echo "[Warn] Missing: ${g}.txt"; continue; }
    awk '{print $2}' "${g}.txt" > "${g}.col2"
    col_files+=("${g}.col2")
done

paste "${col_files[@]}" > gene_all

# Step 4: column-wise max-abs normalization to [-1, 1]
awk '
{
    for (i = 1; i <= NF; i++) {
        data[NR, i] = $i
        absval = ($i < 0) ? -$i : $i
        if (NR == 1 || absval > max[i]) max[i] = absval
    }
    rows = NR
}
END {
    for (r = 1; r <= rows; r++) {
        for (c = 1; c <= NF; c++)
            printf "%f%s", (max[c] > 0 ? data[r,c]/max[c] : 0), (c == NF ? ORS : OFS)
    }
}' gene_all > gene_all_scaled

# Step 5: add gene name header
header=$(awk '{printf "\"%s\",",$3}' "${GENE_LIST}" | sed 's/,$//')
echo "${header}" > gene_t_with_header.txt
cat gene_t_with_header.txt gene_all_scaled > gene_all_scaled_with_header.txt

# Step 6: add cell type row names
printf '"Celltype"\n"BNav"\n"BMem"\n"CD4Naive"\n"CD4EM"\n"CD4Treg"\n"CD8Naive"\n"CD8GZMH"\n"CD8GZMK"\n"MAIT"\n"NKDim"\n"NKBright"\n"MonocM"\n"NoClaM"\n' \
    > first_col.txt

paste -d, first_col.txt gene_all_scaled_with_header.txt \
    | tr ' ' '\t' | sed 's/"//g' | tr ',' '\t' \
    > gene_all_scaled_with_header_with_firstcol.txt

# Step 7: transpose matrix (rows=genes, cols=cell types for pheatmap)
awk '
{
    for (i = 1; i <= NF; i++) a[NR, i] = $i
}
NF > p { p = NF }
END {
    for (j = 1; j <= p; j++) {
        str = a[1, j]
        for (i = 2; i <= NR; i++) str = str "\t" a[i, j]
        print str
    }
}' gene_all_scaled_with_header_with_firstcol.txt > transposed_output.txt

mv transposed_output.txt gene_all_scaled_with_header_with_firstcol.txt

# Step 8: draw heatmap with significance stars
export GENE_LIST
Rscript - << 'REOF'
rm(list = ls()); gc()
library(pheatmap)

infile    <- "gene_all_scaled_with_header_with_firstcol.txt"
outfile   <- "scaled_expression_heatmap_with_stars.pdf"
gene_list <- read.table(Sys.getenv("GENE_LIST"), stringsAsFactors = FALSE)
genes     <- gene_list$V3

cell_order <- c(
  "BNav", "BMem", "CD4Naive", "CD4EM", "CD4Treg",
  "CD8Naive", "CD8GZMH", "CD8GZMK", "MAIT",
  "NKDim", "NKBright", "MonocM", "NoClaM"
)

df <- read.table(infile, header = TRUE, sep = "",
                 check.names = FALSE, stringsAsFactors = FALSE)
rownames(df) <- df[[1]]
df  <- df[, setdiff(colnames(df), colnames(df)[1]), drop = FALSE]
mat <- as.matrix(df[genes, cell_order, drop = FALSE])
mat[mat == 0] <- NA

marker_list <- list(
  c(1,1,1,1,1,1,1,1,1,1,1,1,1),  # HLA-A
  c(1,1,1,1,1,1,1,1,0,1,0,1,1),  # HLA-C
  c(1,1,1,1,0,1,1,1,0,1,1,1,1),  # HLA-F
  c(1,1,0,1,1,0,1,1,0,1,1,1,1),  # HLA-DRB1
  rep(1,13),                       # HLA-DMA
  rep(1,13),                       # GAS5
  c(1,1,1,1,0,1,1,1,0,1,0,1,1),  # ZFAS1
  c(1,1,1,1,1,1,1,1,1,1,0,0,0),  # PCED1B-AS1
  rep(1,13),                       # SNHG15
  rep(1,13),                       # SNHG8
  c(1,1,1,1,1,1,1,1,0,1,0,0,0),  # UQCRB
  c(0,1,1,1,0,1,0,0,0,0,0,0,0),  # ATP5PO
  c(0,1,0,0,0,0,0,0,0,0,0,1,1),  # PRELID1
  c(0,0,0,0,0,0,0,0,0,1,0,0,0),  # TECR
  c(0,1,0,0,0,0,0,1,0,0,0,0,0),  # TYMP
  c(0,0,0,0,0,0,0,0,0,0,0,1,0),  # UBXN1
  c(1,1,1,1,1,1,1,1,1,1,1,0,0),  # NCR3
  c(0,0,1,0,0,0,0,0,0,0,0,0,0),  # IL32
  c(0,0,1,1,0,0,0,0,0,0,0,0,0),  # IL7R
  c(0,0,0,0,0,0,1,1,0,1,1,0,0),  # GNLY
  c(1,0,1,1,0,1,1,1,0,1,0,1,0),  # ITGB2
  c(0,0,0,0,0,0,0,0,0,1,0,0,0),  # CLEC2B
  c(1,1,1,1,0,0,0,0,0,0,0,1,0),  # CD55
  c(1,0,1,1,0,0,0,0,0,0,0,0,0),  # CD37
  c(1,1,1,1,1,0,1,1,0,1,0,1,1),  # RAC1
  c(0,1,1,1,1,0,1,1,0,1,0,1,0),  # CFLAR
  c(0,0,0,0,0,0,0,0,0,0,0,0,1),  # LY96
  c(0,0,1,0,0,0,0,0,0,0,0,0,0),  # RASGRP2
  c(0,0,1,1,0,0,0,0,0,0,0,0,0),  # SIVA1
  c(1,1,1,1,0,1,1,1,0,1,0,0,0),  # FDPS
  c(0,0,1,0,0,0,0,0,0,0,0,0,0),  # SRSF2
  c(1,0,0,0,0,0,0,0,0,0,0,0,0),  # SRSF10
  c(0,0,0,0,0,0,0,0,0,0,0,1,0),  # HNRNPC
  c(0,0,1,1,1,0,1,1,0,1,0,0,0),  # XRN2
  c(1,1,1,1,1,1,1,1,0,1,0,1,1),  # U2AF1L4
  c(0,0,0,0,1,0,0,0,0,0,0,0,0),  # ILF3
  c(1,1,0,0,0,0,0,0,0,1,0,0,0),  # BCL2A1
  c(0,1,0,1,1,0,0,0,0,0,1,0,0),  # CDC42
  c(0,0,1,1,0,1,0,0,0,0,0,0,1)   # NUCKS1
)

data_frame <- do.call(rbind, marker_list)
rownames(data_frame) <- genes
colnames(data_frame) <- cell_order
stars <- ifelse(as.matrix(data_frame[rownames(mat), colnames(mat)]) == 1, "*", "")

pdf(outfile, width = 10, height = max(6, 0.35 * nrow(mat)))
pheatmap(
  mat,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  scale           = "none",
  color           = colorRampPalette(c("blue", "white", "red"))(100),
  breaks          = seq(-1, 1, length.out = 101),
  cellwidth       = 14,
  cellheight      = 14,
  border_color    = "grey85",
  fontsize_row    = 9,
  fontsize_col    = 10,
  display_numbers = stars,
  number_color    = "black",
  na_col          = "grey70"
)
dev.off()
REOF

echo "Done: ${OUT}/scaled_expression_heatmap_with_stars.pdf"