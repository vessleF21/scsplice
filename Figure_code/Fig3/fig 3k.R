#!/usr/bin/env bash

B="/test/3"
FDR="$B/FDR"
NOM="$B/3t/nominal2"
INTRON_OK="$B/3s/intron_13_ok"
N=13

mkdir -p "$B"/{1s,2s/intron,3s,4s/all}

# Step 1: collect intron list per chromosome
for i in {1..22}; do
    cat "$FDR"/*_"$i"_nominals_1.txt \
        | awk '{print $1}' | sort | uniq -c \
        | grep chr | awk '{print $2}' \
        > "$B/1s/intron_chr_$i"
done

# Step 2: count cell types per intron
for i in {1..22}; do
    f="$B/1s/intron_chr_$i"
    [ -f "$f" ] || continue
    while read -r iv; do
        grep "$iv" "$FDR"/*_"$i"_nominals_1.txt \
            | awk -F'/' '{print $NF}' \
            | awk -F'_' '{print $1}' \
            | sort | uniq \
            | tee >(wc -l >> "$B/2s/intron_chr_count_$i") \
            | tr '\n' ',' | sed 's/,$/\n/' \
            >> "$B/2s/intron_chr_cellty_$i"
        grep "$iv" "$FDR"/*_"$i"_nominals_1.txt \
            | awk -F'/' '{print $NF}' \
            | awk -F'[_:]' '{print $1, $5":"$6":"$7"_"$8}' \
            | sort | uniq \
            > "$B/2s/intron/$iv"
    done < "$f"
done

# Step 3: extract shared SNPs across all N cell types
while read -r iv; do
    [ -z "$iv" ] && continue
    c=$(echo "$iv" | awk -F: '{print $1}' | sed 's/chr//')
    grep "$iv" "$NOM"/*_"$c"_nominals_1.txt \
        | awk '{print $11}' | sort | uniq -c \
        | awk -v n=$N '$1==n{print $2}' \
        > "$B/3s/${iv}@13SNP"
done < "$INTRON_OK"

# Step 4: submit PBS jobs
cd "$B/3s/pbs_scripts"
for p in *.pbs; do qsub "$p"; done

# Step 5: aggregate results per intron directory
find "$B/3s" -type d -name "chr*" \
    | awk -F/ '{print $NF}' \
    | while read -r d; do
        base="${d/@13SNP_dir/}"
        mkdir -p "$B/4s/$base"
        cat "$B/3s/$d/"*.txt | awk 'NR%2==1' | sort | uniq -c \
            > "$B/4s/$base/${base}@sign_size"
        cat "$B/3s/$d/"*.txt | awk 'NR%2==0' | sort | uniq -c \
            > "$B/4s/$base/${base}@effect_size"
    done

# Step 6: merge
cd "$B/4s/all"

cat "$B/4s/"*/*effect* \
    | awk '{s[$2]+=$1} END{for(k in s) print s[k],k}' \
    | sort -k2n > effect_2

cat "$B/4s/"*/*sign* \
    | awk '{s[$2]+=$1} END{for(k in s) print s[k],k}' \
    | sort -k2n > sign_2

# Step 7: plot
Rscript - effect_2 sign_2 all.pdf << 'REOF'
args        <- commandArgs(trailingOnly = TRUE)
effect_file <- args[1]
sign_file   <- args[2]
output_file <- args[3]

data_effect <- read.table(effect_file, header = FALSE)
colnames(data_effect) <- c("Count_effect", "Value")
data_effect <- data_effect[order(data_effect$Value), ]
data_effect$percentages_effect <- data_effect$Count_effect / sum(data_effect$Count_effect) * 100

data_sign <- read.table(sign_file, header = FALSE)
colnames(data_sign) <- c("Count_sign", "Value")
data_sign <- data_sign[order(data_sign$Value), ]
data_sign$percentages_sign <- data_sign$Count_sign / sum(data_sign$Count_sign) * 100

merged <- merge(
  data_sign[,  c("Value", "percentages_sign")],
  data_effect[, c("Value", "percentages_effect")],
  by = "Value", all = TRUE
)
merged[is.na(merged)] <- 0

pdf(file = output_file, width = 8, height = 6)

barplot(
  merged$percentages_sign,
  names.arg = merged$Value,
  col       = "orange",
  ylim      = c(-20, 20),
  main      = "Percentage Distribution (Up and Down)",
  xlab      = "Value",
  ylab      = "Percentage (%)",
  border    = "black"
)

barplot(
  -merged$percentages_effect,
  col    = "skyblue",
  add    = TRUE,
  border = "black"
)

dev.off()
REOF