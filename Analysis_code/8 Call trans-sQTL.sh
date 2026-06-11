

##Trans-sQTL Mapping with QTLtools##

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
                                             # Same as cis-sQTL analysis
    \
    --cov PC.txt \                          # Covariate matrix
    \
    --out trans.sample \                    # Output prefix for permutation results
    \
    --sample 5000 \                         # Number of permutations per phenotype
    \
    --threshold 1e-5 \                      # Nominal P-value threshold for reporting
                                             # Only associations P < 1e-5 are saved
                                             # Controls output file size
    \
    --normal \                              # Assume normal distribution of phenotypes
                                             # Valid after inverse normal transformation
    \
    --window 5000000 \                      # Trans-window definition: 5 Mb
                                             # Excludes variants within 5 Mb of phenotype
                                             # Ensures analysis includes only distant effects
    \
    --bin 1000                              # Bin size for computational efficiency

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
# Input: trans.adjust.txt
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
    \
    output                                  # Output prefix for results



# =============================================================================
# Fine-mapping and mediation analysis
# susieR : identifies independent causal signals within a credible set
# mediation: tests whether gene expression mediates the SNP -> splicing effect
#
# Official links:
#   susieR    : https://github.com/stephenslab/susieR
#   mediation : https://cran.r-project.org/package=mediation
# =============================================================================
library(susieR)
library(mediation)

# =============================================================================
# Part 1. Fine-mapping with SuSiE RSS (summary statistics input)
# Uses z-scores and an LD matrix to compute posterior inclusion probabilities
# (PIPs) for each variant and identify credible sets of causal variants.
#
# Input:
#   z  : z-scores computed as beta / se from cis-eQTL or cis-sQTL results
#   R  : in-sample LD matrix (correlation matrix from the same genotype data)
#   n  : sample size
#   L  : maximum number of independent causal signals to consider
#
# Output:
#   result$pip      : posterior inclusion probability per variant
#   result$sets     : credible sets at specified coverage
# =============================================================================
result <- susie_rss(
  z        = data$zscore,
  R        = R,
  n        = nrow(R),
  L        = 2,
  estimate_residual_variance = TRUE,
  coverage = 0.95
)

# =============================================================================
# Part 2. Mediation analysis
# Tests whether gene expression (M) mediates the effect of a genetic variant
# (X) on splicing (Y). Decomposes total effect into:
#   ACME  : average causal mediation effect (indirect: X -> M -> Y)
#   ADE   : average direct effect           (direct:   X -> Y)
#   Total : total effect = ACME + ADE
#   Prop. Mediated: proportion of total effect explained by mediation
#
# Input:
#   X : genotype dosage (0/1/2)
#   M : gene expression (mediator)
#   Y : PSI / splicing ratio (outcome)
#
# Bootstrap (sims=1000) used to estimate confidence intervals for ACME.
# =============================================================================
mdl_m <- lm(M ~ X, data = df)          # alpha path: X -> M
mdl_y <- lm(Y ~ X + M, data = df)      # beta path:  M -> Y (adjusting for X)

med_out <- mediate(mdl_m, mdl_y,
                   treat    = "X",
                   mediator = "M",
                   boot     = TRUE,
                   sims     = 1000)
summary(med_out)







    