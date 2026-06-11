
##This pipeline implements three complementary QTLtools modes to comprehensively identify genetic variants affecting splicing##

# Official qtltools repository link: https://qtltools.github.io/qtltools/

# Permutation Mode: Identify top cis-sQTL per gene (corrected P-values)
# Conditional Mode: Discover independent cis-sQTL signals
# Nominal Mode: Detect all variant-junction associations


# Step 1: Permutation Pass - Top cis-sQTL Identification
##=============================================================================
## PERMUTATION PASS MODE
##=============================================================================
# Identify the top nominal cis-sQTL for each cis-sGene with empirical P-value correction
#
# Purpose:
# - Establishes phenotype-specific null distribution via permutations
# - Corrects for multiple testing across variants within cis-window
# - Identifies most significant variant per junction cluster/gene

QTLtools cis \
    --vcf filter.chr.dose.vcf.gz \              # Imputed genotypes (dosage format)

    \
    --bed celltype_chr_phenotype.txt.gz \       # Normalized junction usage (phenotypes)

    \
    --cov celltype_ok_PC.txt \                  # Covariate matrix
    \
    --permute 1000 \                            # Number of permutations for empirical P-value
                                            
    \
    --normal \                                  # Assume normal distribution for phenotypes
                                                 # Appropriate for normalized PSI values
    \
    --grp-best \                                # Group junctions by gene (cluster-level analysis)
                                                 # Aggregates associations across related junctions
                                                 # Calculates empirical P-value per gene
                                                 # Implements: "comprehensive assessment across 
                                                 # junction clusters linked to a given gene"
    \
    --out celltype_chr_permutation.txt          # Output

# Run for each chromosome separately (parallelization)

##-----------------------------------------------------------------------------
## Aggregate Results Across Chromosomes
##-----------------------------------------------------------------------------
# Combine chromosome-specific permutation results into single file
cat celltype_*_permutation.txt | gzip -c > celltype_permutations_full.txt.gz

##-----------------------------------------------------------------------------
## FDR Correction
##-----------------------------------------------------------------------------
# Apply FDR correction
# Implements: "Statistical significance determined using FDR threshold"

Rscript runFDR_cis.R \
    celltype_permutations_full.txt.gz \         # Input: aggregated permutation results
    0.05 \                                      # FDR threshold (5%)
    sig_celltype_permutations_full.txt          # Output: significant sQTLs only

# Output also includes:
#   - *.thresholds.txt: Phenotype-specific nominal P-value thresholds
#   - Used as input for conditional analysis

#Step 2: Conditional Pass - Independent Signal Discovery
##=============================================================================
## CONDITIONAL PASS MODE
##=============================================================================
# Identify independent cis-sQTL signals per gene using stepwise regression
#
# Purpose:
# - Detect secondary signals after conditioning on lead variant
# - Account for LD structure between variants
# - Provide complete picture of regulatory architecture per gene
#
# Algorithm:
# 1. Start with top variant from permutation pass
# 2. Condition on this variant
# 3. Test remaining variants for residual associations
# 4. If significant, add as independent signal
# 5. Repeat until no more significant signals
# 6. Backward pass: verify all signals remain significant

QTLtools cis \
    --vcf filter.chr.dose.vcf.gz \              # Same genotype input
    \
    --bed celltype_chr_phenotype.txt.gz \       # Same phenotype input
    \
    --cov celltype_ok_PC.txt \                  # Same covariates
    \
    --normal \                                  # Normal distribution assumption
    \
    --grp-best \                                # Gene-level aggregation
    \
    --mapping celltype_permutations_full.txt.thresholds.txt \
                                                 # Phenotype-specific P-value thresholds
    \
    --out conditions_celltype_chr.txt           # Output

# Output format includes:
# - Primary signal (from permutation pass)
# - Secondary signals (after conditioning on primary)
# - Tertiary signals (after conditioning on primary + secondary)
# - etc.
#
# Each line represents an independent association

#Step 3: Nominal Pass - Complete Association Catalog
##=============================================================================
## NOMINAL PASS MODE
##=============================================================================
# Detect all variant-junction associations within cis-window
#
# Purpose:
# - Comprehensive catalog of all associations
# - No multiple testing correction (reports all P-values)

QTLtools cis \
    --vcf filter.chr.dose.vcf.gz \              # Genotype input
    \
    --bed celltype_chr_phenotype.txt.gz \       # Phenotype input
    \
    --cov celltype_ok_PC.txt \                  # Covariate matrix
    \
    --nominal 1 \                               # Report all associations (no filtering)
                                                 # 1 = output all variant-phenotype pairs
                                                 # 0.01 = output pairs with P < 0.01
                                                 # Implements: "detect all available cis-sQTLs"
    \
    --normal \                                  # Normal distribution
    \
    --grp-best \                                # Gene-level grouping
    \
    --out celltype_chr_nominals_1.txt           # Output: complete association table



# =============================================================================
# Part D. RBP Motif Analysis for cis-sQTL Mechanism
# ±50 bp sequences flanking independent cis-sQTL variants (GENCODE-annotated
# introns only) were scanned using FIMO (v5.5.8) with ATtRACT PWMs.
# REF and ALT allele sequences were compared to identify allele-specific RBP
# binding changes. Hypergeometric tests identified cell-type-specific
# RBPs, whose representative motifs were constructed from binding site sequences.
#
# ATtRACT : https://attract.cnic.es
# FIMO    : https://meme-suite.org/meme/tools/fimo
# bedtools: https://github.com/arq5x/bedtools
# =============================================================================

CELLS=(BMem BNav CD4EM CD4Naive CD4Treg CD8GZMH CD8GZMK CD8Naive MAIT MonocM NKBright NKDim NoClaM)

awk '{print $1":"$2":"$3":"$4}' gencode32_introns.txt > pos_strand.txt

# D1. Extract ±50 bp sequences flanking each sQTL SNP from hg38
# snp.bed is generated from GENCODE-filtered sQTL SNP positions (±50 bp window)
# Output: snp_ref.fa contains the reference allele sequences for each SNP
# Official bedtools repository: https://github.com/arq5x/bedtools
for cell in "${CELLS[@]}"; do
    cd $cell
    awk -v ct="$cell" '$23==ct' qtl.txt | \
        awk 'FNR==NR{keep[$0];next}($1 in keep)' pos_strand.txt - > qtl_f.txt
    awk '{print $11}' qtl_f.txt | \
        awk -F: '{s=$2-50;if(s<0)s=0;print $1"\t"s"\t"$2+50"\t"$0}' > snp.bed
    bedtools getfasta \
        -fi hg38.fa \
        -bed snp.bed \
        -fo snp_ref.fa \
        -nameOnly

# D2. Generate REF and ALT allele FASTA files
# allele_seq_generator.py replaces the central base of each sequence with REF or ALT allele
# For minus-strand introns, sequences are reverse-complemented
# Input : snp_ref.fa (header format: chr:pos:ref:alt:strand)
# Output: one *_ref and one *_alt FASTA file per SNP
cat > allele_seq_generator.py << 'PYEOF'
def rc(s): return s.translate(str.maketrans("ACGTacgt","TGCAtgca"))[::-1]
def run(fa):
    lines = open(fa).readlines()
    for i in range(0, len(lines), 2):
        h = lines[i].strip(); s = lines[i+1].strip()
        if not h.startswith('>'): continue
        p = h[1:].split(':')
        if len(p) < 5: continue
        _, _, ref, alt, strand = p[:5]; sid = h[1:]; mid = (len(s)-1)//2
        if s[mid].upper() not in (ref, alt): continue
        for al, b in [('ref',ref),('alt',alt)]:
            seq = (rc(s[:mid]+b+s[mid+1:]) if strand=='-' else s[:mid]+b+s[mid+1:]).upper()
            open(f"{sid}_{al}",'w').write(f">{sid}_{al}\n{seq}\n")
run("snp_ref.fa")
PYEOF

python allele_seq_generator.py


# D3. FIMO motif scanning on REF and ALT sequences
# Scans each allele-specific FASTA against ATtRACT RBP PWMs
# --thresh 0.001: p-value threshold for motif hit reporting
# Output: fimo.tsv containing all significant RBP-sequence matches per SNP
# Official FIMO documentation: https://meme-suite.org/meme/tools/fimo
    for allele in ref alt; do
        for fa in ${allele}_fa/chr*; do
            fimo \
                --thresh 0.001 \
                --oc ${allele}_fa/fimo_results/$(basename $fa) \
                ATtRACT_Homo_sapiens.meme \
                $fa
        done
    done

# D4. Compare REF vs ALT motif hits to identify allele-specific RBP binding
# ref_only: RBP motif present in REF but lost in ALT (SNP disrupts binding)
# alt_only: RBP motif absent in REF but gained in ALT (SNP creates new binding)
    for af in alt_fa/fimo_results/*/fimo_parsed.tsv; do
        rf=${af/alt_fa/ref_fa}; rf=${rf/_alt/_ref}
        [[ -f $rf ]] || continue
        out=result_fa/$(basename $(dirname $af)).tsv; mkdir -p result_fa
        comm -23 <(awk '{print $2}' $rf|sort -u) <(awk '{print $2}' $af|sort -u) | \
            awk '{print $0"\tref_only"}' >  $out
        comm -13 <(awk '{print $2}' $rf|sort -u) <(awk '{print $2}' $af|sort -u) | \
            awk '{print $0"\talt_only"}' >> $out
    done

    cd 
done

# D5. Hypergeometric test for cell-type-specific RBP enrichment (R)
# Hypergeometric test comparing RBP hit frequency in one cell type
# vs all others; RBPs also in top-3 of other cell types are excluded.
# Hit counts are normalized by sQTL number and log10-transformed for heatmap.

Rscript - << 'RSCRIPT'
library(pheatmap)

# One-sided hypergeometric test: x = RBP hits in focal cell, m = total RBP hits,
# n = total hits not in focal cell, k = total hits in focal cell
p_hyper <- phyper(x-1, m=total_rbp, n=M-total_rbp, k=cell_total, lower.tail=FALSE)
p_fisher <- fisher.test(matrix(c(x, cell_total-x, m, other_total-m), 2),
                        alternative="greater")$p.value

# Normalize hit counts by number of sQTLs per cell type; log10 transform
normalized <- log10(rbp_count / n_sqtl)

# Plot heatmap grouped by RBP functional category (white = low, red = high)
pheatmap(mat,
         cluster_rows   = FALSE,
         cluster_cols   = FALSE,
         color          = colorRampPalette(c("white","red"))(100),
         breaks         = seq(0, 1, length.out=101),
         na_col         = "black",
         cellheight     = 10,
         cellwidth      = 10,
         border_color   = "black")
RSCRIPT


# =============================================================================
# Part E. Sex-Biased cis-sQTL Analysis
# Model: Y = B0 + Bg·G + Bd·D + Bgxd·GxD + Ba·A +
#            sum(BsPC·splicingPC) + sum(BgPC·GenotypePC)
#   Y: PSI  G: genotype(0/1/2)  D: sex  A: age
#   splicingPC1-7 / GenotypePC1-5
# Multiple testing:
#   Stage 1: Bonferroni across SNPs per junction
#   Stage 2: Multiple testing (threshold 0.25)
# =============================================================================

# Load qvalue package for FDR estimation
library(qvalue)

# Define 13 immune cell types to analyze
cell_types <- c("CD4EM","CD4Naive","CD4Treg","CD8GZMH","CD8GZMK","CD8Naive",
                "MonocM","MAIT","BMem","BNav","NoClaM","NKBright","NKDim")

for (ctype in cell_types) {

  # Extract covariate matrix: rows 1-12,14 = splicingPC1-7 + GenotypePC1-5 + age
  # Row 13 = sex (D), extracted separately as the interaction term
  cov <- setNames(as.data.frame(t(covar_mx[c(1:12,14),])),
                  c(paste0("spl",1:7), paste0("gtp",1:5), "A"))
  D   <- as.numeric(t(covar_mx[13,]))

  # Fit interaction model per junction-SNP pair
  # Y ~ G*D expands to Y ~ G + D + G:D; G:D tests sex-biased sQTL effect
  # Extract p-value for the G:D interaction term
  res <- na.omit(do.call(rbind, lapply(seq_len(nrow(top_qtl)), function(r) {
    cf <- summary(lm(Y ~ G*D + spl1+spl2+spl3+spl4+spl5+spl6+spl7+gtp1+gtp2+gtp3+gtp4+gtp5+A,
                     data=data.frame(Y=as.numeric(t(psi_mx[r,])),
                                     G=as.numeric(t(geno_mx[r,])), D=D, cov)))$coefficients
    data.frame(junction=top_qtl[r,6], snp=top_qtl[r,10], p=cf["G:D","Pr(>|t|)"])
  })))

  # Stage 1: Bonferroni correction across SNPs within each junction
  res$p <- ave(res$p, res$junction, FUN=function(p) p.adjust(p, method="bonferroni"))

  # Stage 2: For each junction, retain the SNP with the smallest Bonferroni p-value
  best  <- do.call(rbind, lapply(split(res, res$junction), function(s) s[which.min(s$p),]))

  # Stage 2: Apply Storey q-value FDR across junctions; retain significant junctions (q <= 0.25)
  sig   <- best[qvalue(best$p, fdr.level=0.25, pi0=1)$qvalues <= 0.25, "junction"]

  # Write all SNP-junction pairs for significant junctions to output file
  write.table(res[res$junction %in% sig,], paste0("sig_pairs_",ctype,".txt"),
              sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}





