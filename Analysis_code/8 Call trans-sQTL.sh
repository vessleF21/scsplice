

Trans-sQTL Mapping with QTLtools

# Official qtltools repository link: https://qtltools.github.io/qtltools/

#Step 1: Permutation-based Beta Distribution Approximation
##=============================================================================
## TRANS-sQTL MAPPING - PERMUTATION PASS
##=============================================================================
# Calculate empirical P-values through permutation testing
# Establishes beta distribution parameters for efficient trans-sQTL discovery
#

QTLtools trans \
    --vcf chr1_22.0.9.vcf.gz \              # Filtered genotype file
                                             # Additional filtering applied:
                                             # - 50-mer mappability >= 0.9
                                             # - Excludes poorly mappable regions
                                             # - Reduces mapping artifacts
                                             # - Based on k50.umap.bedgraph.gz
                                             # 
    \
    --bed phenotype.txt.gz \                # Normalized splicing matrix
                                             # Same as cis-sQTL analysis:
                                             # - Junction usage (PSI values)
    \
    --cov PC.txt \                          # Covariate matrix
                                             # Same as cis-sQTL analysis:
                                             # - 7 splicing PCs
                                             # - 5 genetic PCs
                                             # - Sex
                                             # - Age
    \
    --out trans.sample \                    # Output prefix for permutation results
    \
    --sample 5000 \                         # Number of permutations per phenotype
                                             # Creates empirical null distribution
                                             # Higher N = more accurate beta parameters
    \
    --threshold 1e-5 \                      # Nominal P-value threshold for reporting
                                             # Only associations P < 1e-5 are saved
                                             # Controls output file size
                                             # Stringent to manage multiple testing
    \
    --normal \                              # Assume normal distribution of phenotypes
                                             # Valid after inverse normal transformation
    \
    --window 5000000 \                      # Trans-window definition: 5 Mb
                                             # Excludes variants within 5 Mb of phenotype
                                             # Ensures analysis includes only distant effects
    \
    --bin 1000                              # Bin size for computational efficiency
                                             # Groups nearby variants for faster processing
                                             # Does not affect statistical validity

# Output files:
# - trans.sample.txt: Summary statistics for each phenotype
#   Columns:
#   - phenotype_id: Junction cluster ID
#   - n_variants: Number of variants tested
#   - n_permutations: Number of permutations performed
#   - beta_shape1: First beta distribution parameter (α)
#   - beta_shape2: Second beta distribution parameter (β)
#   - true_df: Degrees of freedom
#

##=============================================================================
## PERMUTATION DETAILS
##=============================================================================
# For each phenotype:
# 1. Test all trans variants against observed phenotype
# 2. Repeat 5000 times with permuted sample labels
# 3. Record minimum P-value from each permutation
# 4. Fit beta distribution to these minimum P-values
# 5. Beta parameters (α, β) characterize the null distribution
#
# Beta approximation advantages:
# - Fast: No need to re-permute in adjustment step
# - Accurate: Captures phenotype-specific null
# - Scalable: Enables genome-wide trans analysis

#Step 2: Multiple Testing Adjustment
##=============================================================================
## TRANS-sQTL MAPPING - ADJUSTMENT PASS
##=============================================================================
# Adjust P-values for multiple testing across phenotypes per gene


QTLtools trans \
    --vcf chr1_22.0.9.vcf.gz \              # Same genotype file (mappability filtered)
    \
    --bed phenotype.txt.gz \                # Same phenotype file
    \
    --cov PC.txt \                          # Same covariates
    \
    --out trans.adjust \                    # Output prefix for adjusted results
    \
    --adjust trans.best.txt.gz \            # Input: Best trans associations
    \
    --threshold 0.1 \                       # Adjusted P-value threshold
    \
    --normal \                              # Normal distribution assumption
    \
    --window 5000000 \                      # Same 5 Mb trans-window
                                             # Consistent with permutation step
    \
    --bin 1000                              # Same bin size

# Output files:
# - trans.adjust.txt: Adjusted P-values for trans associations
#   Columns:
#   - phenotype_id: Junction cluster ID
#   - variant_id: Trans-acting variant
#   - chr_variant: Variant chromosome
#   - pos_variant: Variant position
#   - chr_phenotype: Phenotype chromosome
#   - pos_phenotype: Phenotype position
#   - distance: Genomic distance (>5Mb or cross-chromosome)
#   - pval_nominal: Original nominal P-value
#   - pval_adjusted: Empirically adjusted P-value
#   - slope: Effect size (beta coefficient)
#   - slope_se: Standard error
#


Step 3: FDR Correction Across All Trans-sQTLs
##=============================================================================
## GENOME-WIDE FDR CORRECTION
##=============================================================================
# Apply FDR correction across all trans-sQTL candidates
#
# Input: trans.adjust.txt (from adjustment step)
# Process: FDR correction considering mappability
# Output: Significant trans-sQTLs (FDR < 0.05)

Rscript runFDR.mappability.R \
    celltype_chr_trans.txt.gz \             # Input: Adjusted trans-sQTL results
    \
    0.05 \                                  # FDR threshold (5%)
                                             # Standard threshold for QTL studies
                                             # Controls false discovery rate
    \
    trans.adjust.txt \                      # Adjustment statistics file
                                             # Contains beta parameters and correction info
    \
    output                                  # Output prefix for results