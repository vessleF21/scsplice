##This pipeline implements three complementary QTLtools modes to comprehensively identify genetic variants affecting splicing##

# Official qtltools repository link: https://qtltools.github.io/qtltools/

# Permutation Mode: Identify top cis-sQTL per gene (corrected P-values)
# Conditional Mode: Discover independent cis-sQTL signals
# Nominal Mode: Detect all variant-junction associations (complete catalog)


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
                                                 # Each phenotype permuted 1000 times to build null
                                                 # Higher N = more accurate empirical P-values
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
    --out celltype_chr_permutation.txt          # Output: one line per phenotype with top variant

# Run for each chromosome separately (parallelization)
# Output format:
# phenotype_id | chr | variant_pos | variant_id | ... | beta | empirical_pval | adj_beta_pval

##-----------------------------------------------------------------------------
## Aggregate Results Across Chromosomes
##-----------------------------------------------------------------------------
# Combine chromosome-specific permutation results into single file
cat celltype_*_permutation.txt | gzip -c > celltype_permutations_full.txt.gz

##-----------------------------------------------------------------------------
## FDR Correction
##-----------------------------------------------------------------------------
# Apply Storey-Tibshirani FDR correction across all phenotypes
# Implements: "Statistical significance determined using FDR threshold"
# 
# This script:
# 1. Reads empirical P-values from permutation pass
# 2. Applies Benjamini-Hochberg or q-value FDR correction
# 3. Determines phenotype-specific P-value thresholds at FDR = 0.05
# 4. Outputs significant cis-sGenes (FDR < 0.05)

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
# 2. Condition on this variant (include as covariate)
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
                                                 # From FDR correction step
                                                 # Determines significance for each phenotype
                                                 # Implements: "phenotype-specific nominal P value 
                                                 # threshold, which varied according to the number 
                                                 # of independent tests"
    \
    --out conditions_celltype_chr.txt           # Output: all independent signals per phenotype

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