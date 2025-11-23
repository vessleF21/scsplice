

##TensorQTL-based Cis-eQTL Mapping Pipeline##

# Official tensorQTL repository link: https://github.com/broadinstitute/tensorqtl

#Import Dependencies and Setup
##=============================================================================
## IMPORT LIBRARIES AND INITIALIZE
##=============================================================================
import sys
import rpy2.robjects as ro
import tensorqtl
import pandas as pd
import torch
import rpy2
from tensorqtl import pgen, cis, trans, post
from pyplink import PyPlink

# Parse command-line argument: cell type identifier
# Enables parallel processing across cell types
celltype = sys.argv[1]

# Define input file paths for eQTL analysis
# Note: Different from sQTL - using gene expression instead of junction usage
phenotype_bed_file = f"{celltype}_sorted.bed.gz"     # Pseudobulk gene expression
covariates_file = f"{celltype}_ok_PC.txt"             # Covariate matrix (with PEER factors)
pgen_file = f"{celltype}_tsst"                        # PLINK2 genotype format

#Step 1: Load Input Data
##=============================================================================
## LOAD PHENOTYPE DATA - PSEUDOBULK GENE EXPRESSION
##=============================================================================
# Read normalized gene expression levels (BED format)
#
# Phenotype preparation:
# 1. Generation: "pseudobulk expression matrices by averaging UMI counts 
#    across all cells of each cell type per sample"
# 2. Filtering: "Genes detected in fewer than 40% of individuals were excluded"
# 3. Normalization:
#    a. "quantile normalization with trimmed mean of M-values (TMM) adjustment"
#       - Normalizes across samples (removes library size effects)
#    b. "inverse quantile transformation to a standard normal distribution"
#       - Normalizes within genes (rank-based normalization)
#       - Ensures normal distribution assumption for linear regression
#
# BED format:
# - chr, start, end, gene_id, sample1, sample2, ...
# - Values: Normalized gene expression (standard normal distribution)
# - TSS (transcription start site) used as gene position

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)

# phenotype_df: Sample × Gene matrix (normalized expression values)
# phenotype_pos_df: Gene coordinates (chr, start=TSS, end=TES, gene_id)

##=============================================================================
## LOAD COVARIATE DATA - INCLUDING PEER FACTORS
##=============================================================================
# Read covariate matrix and transpose to Sample × Covariate format
#
# Implements: "Covariates included age, sex, 5 genotype PCs, and probabilistic 
# estimation of expression residuals (PEER) factors"
#
#
# PEER factors (Probabilistic Estimation of Expression Residuals):
# - Capture hidden confounders in gene expression data
# - Technical variation: batch effects, library prep, sequencing depth
# - Biological variation: cell composition, cell cycle, stress response
#
# Number of PEER factors by sample size (from methodology):
# "The number of PEER factors per cell type was selected according to 
# sample number: 15 for N < 150, 30 for 150 ≤ N < 250, 45 for 250 ≤ N < 350, 
# and 60 for N ≥ 350"
#

covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

##=============================================================================
## LOAD GENOTYPE DATA
##=============================================================================
# Load genotypes in PLINK2 PGEN format
# 
# Implements: "Individuals with fewer than 10 cells per cell type were excluded"
#
# Genotype QC:
# - Call rate >= 95%
# - MAF >= 5% (after imputation)
# - HWE P-value > 1e-6
# - Imputation quality R² > 0.8
# - Biallelic autosomal SNPs only
#

pgr = pgen.PgenReader(pgen_file)
genotype_df = pgr.load_genotypes()      # Variant × Sample matrix (dosage format)
variant_df = pgr.variant_df              # Variant annotation dataframe

# Standardize chromosome naming to match phenotype file
variant_df['chrom'] = variant_df['chrom'].apply(lambda x: 'chr' + str(x))

#Step 2: Cis-eQTL Mapping (Permutation Mode)
##=============================================================================
## CIS-eQTL MAPPING - PERMUTATION PASS
##=============================================================================
# Identify top cis-eQTL per gene with empirical P-value correction
#
# Purpose: Link genetic variants to gene expression levels
#
# Analysis parameters:
# - Cis-window: 1 Mb
# - Model: Linear regression (genotype dosage ~ normalized expression)
# - Covariates: Genetic PCs, sex, age, PEER factors
# - Permutations: Empirical null distribution per gene
#
# Statistical model:
# Expression_g = β₀ + β₁·Genotype_v + β₂·gPC + β₃·PEER + β₄·age + β₅·sex + ε

cis_df = cis.map_cis(
    genotype_df,         # Variant × Sample genotype matrix
    variant_df,          # Variant annotations (chr, pos, id, ref, alt)
    phenotype_df,        # Sample × Gene matrix (normalized expression)
    phenotype_pos_df,    # Gene coordinates (TSS as primary position)
    covariates_df        # Sample × Covariate matrix (gPC + PEER + age + sex)
)

# Output columns in cis_df:
# - phenotype_id: Gene identifier (ENSG ID or gene symbol)
# - variant_id: Top associated variant (lead eQTL)
# - tss_distance: Distance from variant to gene TSS
# - pval_nominal: Nominal P-value (uncorrected)
# - slope: Effect size (beta coefficient)
#   * Positive: Allele increases expression
#   * Negative: Allele decreases expression
# - slope_se: Standard error of effect size
# - pval_perm: Empirical P-value from permutations
# - pval_beta: Beta-approximated P-value

##=============================================================================
## Q-VALUE CALCULATION
##=============================================================================
# Calculate q-values for multiple testing correction
#
# Implements: "Significant cis-eQTLs were defined at FDR threshold of < 0.05 
# after Benjamini-Hochberg correction"
#
# Two-step FDR correction procedure:

# Step 1: Calculate q-values across all genes
# qvalue_lambda=0.85: Tuning parameter for π₀ estimation
# - π₀ = proportion of true null hypotheses
# - lambda controls robustness of estimation
# - Higher lambda = more conservative
tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0.85)

# Step 2: Apply FDR threshold (5%) and filter significant eQTLs
# Adds 'qval' column to cis_df
# Genes with qval < 0.05 are significant eGenes (expression QTL genes)
post.calculate_qvalues(cis_df, fdr=0.05, qvalue_lambda=0.85)

# Save cis-eQTL results with q-values
# Contains lead eQTL (top variant) for each significant eGene
cis_df.to_csv(f'{celltype}_cis_df_saved.txt', sep='\t', index=True)

#Step 3: Nominal Pass - Complete Association Catalog
##=============================================================================
## CIS-eQTL MAPPING - NOMINAL PASS
##=============================================================================
# Report all variant-gene associations within cis-window
#
# Purpose:
# - Comprehensive catalog for fine-mapping
# Output: Chromosome-specific files for efficient storage

prefix = f"{celltype}_cisqtl_sample"

cis.map_nominal(
    genotype_df,         # Genotype matrix
    variant_df,          # Variant annotations
    phenotype_df,        # Gene expression matrix
    phenotype_pos_df,    # Gene TSS coordinates
    prefix,              # Output file prefix
    covariates_df,       # Covariate matrix (with PEER factors)
    output_dir='.'       # Output directory
)

# Output files: {prefix}.cis_qtl_pairs.chr{1-22}.parquet
# Each file contains all variant-gene pairs for one chromosome
# Columns:
# - phenotype_id: Gene ID
# - variant_id: Variant ID
# - tss_distance: Distance from variant to TSS
# - pval_nominal: Nominal P-value (no correction)
# - slope: Effect size
# - slope_se: Standard error

##=============================================================================
## CONVERT PARQUET TO TEXT FORMAT
##=============================================================================
# Convert binary parquet files to tab-delimited text
# Facilitates downstream colocalization and integration analyses

parquet_files = [f"{prefix}.cis_qtl_pairs.chr{i}.parquet" for i in range(1, 23)]

for parquet_file in parquet_files:
    try:
        # Read parquet file
        data = pd.read_parquet(parquet_file)
        
        # Write as tab-delimited text
        txt_file = parquet_file.replace('.parquet', '.txt')
        data.to_csv(txt_file, sep='\t', index=False)
        
        print(f"Converted {parquet_file} to {txt_file}")
    except Exception as e:
        print(f"Failed to process {parquet_file}: {e}")

#Step 4: Independent Signal Discovery (Conditional Analysis)
##=============================================================================
## INDEPENDENT SIGNAL MAPPING - CONDITIONAL ANALYSIS
##=============================================================================
# Identify independent cis-eQTL signals using stepwise regression
#
# Significance: Bonferroni correction within cis-window
# Accounts for LD between variants and effective number of tests

indep_df = cis.map_independent(
    genotype_df,         # Genotype matrix
    variant_df,          # Variant annotations
    cis_df,              # Results from permutation pass (significant eGenes)
    phenotype_df,        # Gene expression matrix
    phenotype_pos_df,    # Gene TSS coordinates
    covariates_df        # Covariate matrix
)

# Output columns in indep_df:
# - phenotype_id: Gene identifier
# - variant_id: Independent eQTL variant
# - rank: Signal rank (1=primary, 2=secondary, 3=tertiary, etc.)
# - pval_nominal: Nominal P-value after conditioning on previous signals
# - slope: Effect size (independent of other signals)
# - tests_emt: Effective number of independent tests (multiple testing burden)
#

# Save independent signal results
indep_df.to_csv(f'{celltype}_indep_df_saved.txt', sep='\t', index=True)