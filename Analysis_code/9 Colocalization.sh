

##sQTL-GWAS Colocalization Analysis using COLOC##

# Official coloc repository link: https://github.com/chr1swallace/coloc

#Define GWAS Traits and Study Parameters
##=============================================================================
## GWAS TRAIT DEFINITIONS - CASE-CONTROL STUDIES
##=============================================================================
# Define case-control GWAS traits for colocalization analysis
# Format: Author_Journal_Year_Disease
# Covers autoimmune diseases, cancers, and metabolic disorders

cc_trait <- c(
    "ValetteK_CommunBiol_2021_asthma",           # Respiratory
    "LiuZ_NatGenet_2023_CD",                     # Crohn's disease
    "AndlauerTF_SciAdv_2016_MS",                 # Multiple sclerosis
    "Ishigaki_NatGenet_2022_RA",                 # Rheumatoid arthritis
    "BenthamJ_NatGenet_2015_SLE",                # Systemic lupus erythematosus
    "SakaueS_NatGenet_2021_SS",                  # Sjögren's syndrome
    "ChiouJ_Nature_2021_T1D",                    # Type 1 diabetes
    "LiuZ_NatGenet_2023_UC",                     # Ulcerative colitis
    "SakaueS_NatGenet_2021_GD",                  # Graves' disease
    "MichailidouK_Nature_2017_Breastcancer",     # Breast cancer
    "WangA_NatGenet_2023_Prostatecancer",        # Prostate cancer
    "McKayJD_NatGenet_2017_lungcarcinoma",       # Lung cancer
    "PurdueMP_NatGenet_2024_Kidneycancer",       # Kidney cancer
    "SeviiriM_NatCommun_2022_Basalcellcarcinoma", # Basal cell carcinoma
    "MaraTA_NatCommun_2018_endometrialcarcinoma", # Endometrial cancer
    "SeviiriM_NatCommun_2022_squamouscellcarcinoma", # Squamous cell carcinoma
    "RamachandranD_HumMolGenet_2022_cervicalcancer", # Cervical cancer
    "DarengEO_AmJHumGenet_2024_ovariancancer",   # Ovarian cancer
    "SakaueS_NatGenet_2021_Type2diabetes",       # Type 2 diabetes
    "TrinderM_Atherosclerosis_2021_hyperlipidemia" # Hyperlipidemia
)

# Sample sizes for case-control studies
# Critical for accurate colocalization calculation (affects prior probabilities)
cc_trait_case <- c(56167, 20873, 4888, 22350, 5201, 1599, 18942, 23252, 4487, 
                   76192, 122188, 29266, 29020, 20791, 8758, 7402, 363, 15588, 
                   84224, 17485)

cc_trait_control <- c(352255, 346719, 10395, 74823, 9066, 658316, 501638, 
                      352256, 629598, 63082, 604640, 56450, 835670, 286893, 
                      46126, 286892, 861, 105724, 583280, 331737)

##=============================================================================
## GWAS TRAIT DEFINITIONS - QUANTITATIVE TRAITS
##=============================================================================
# Define quantitative GWAS traits
# Includes blood cell counts, anthropometric traits, metabolic phenotypes

quant_trait <- c(
    "ChenMH_CELL_2020_BAS",                      # Basophil count
    "SakaueS_NatGenet_2021_BMI",                 # Body mass index
    "ChenMH_CELL_2020_EOS",                      # Eosinophil count
    "SakaueS_NatGenet_2021_height",              # Height
    "ChenMH_CELL_2020_Ht",                       # Hematocrit
    "ChenMH_CELL_2020_LYM",                      # Lymphocyte count
    "ChenMH_CELL_2020_MCHC",                     # Mean corpuscular hemoglobin concentration
    "ChenMH_CELL_2020_MCH",                      # Mean corpuscular hemoglobin
    "ChenMH_CELL_2020_MCV",                      # Mean corpuscular volume
    "ChenMH_CELL_2020_MON",                      # Monocyte count
    "ChenMH_CELL_2020_NEU",                      # Neutrophil count
    "ChenMH_CELL_2020_PLT",                      # Platelet count
    "ChenMH_CELL_2020_RBC",                      # Red blood cell count
    "ChenMH_CELL_2020_WBC",                      # White blood cell count
    "ChenMH_CELL_2020_Hb",                       # Hemoglobin
    "ParkS_NatGenet_2024_Metabolicsyndrome",     # Metabolic syndrome
    "CareyCE_NatHumBehav_2024_Hypertension"      # Hypertension
)

# Sample sizes for quantitative traits
quant_trait_num <- c(474001, 523818, 474237, 525444, 562259, 524923, 491553, 
                     486823, 544127, 521594, 519288, 542827, 545203, 562243, 
                     563946, 1384348, 338391)

##=============================================================================
## CELL TYPE DEFINITIONS AND SAMPLE SIZES
##=============================================================================
# Define cell types analyzed in sQTL study
# Sample sizes critical for COLOC analysis (affects statistical power)

cell_type <- c("BMem",        # Memory B cells
               "BNav",        # Naive B cells
               "CD4EM",       # CD4+ effector memory T cells
               "CD4Naive",    # CD4+ naive T cells
               "CD4Treg",     # CD4+ regulatory T cells
               "CD8GZMH",     # CD8+ GZMH+ T cells
               "CD8GZMK",     # CD8+ GZMK+ T cells
               "CD8Naive",    # CD8+ naive T cells
               "MAIT",        # Mucosal-associated invariant T cells
               "MonocM",      # Classical monocytes
               "NKBright",    # CD56bright NK cells
               "NKDim",       # CD56dim NK cells
               "NoClaM")      # Non-classical monocytes

# Number of samples per cell type (after QC filtering)
# Varies due to cell type abundance and quality thresholds
celltype_sample <- c(865, 873, 980, 975, 729, 939, 924, 759, 352, 603, 
                     346, 976, 435)

#Parse Command-line Arguments
##=============================================================================
## COMMAND-LINE ARGUMENT PARSING
##=============================================================================
# Script expects 4 arguments for parallel processing across traits/loci
# Enables high-throughput colocalization testing

args <- (commandArgs(TRUE))
gwas_name <- args[1]    # GWAS trait name (e.g., "LiuZ_NatGenet_2023_CD")
celltype <- args[2]     # Cell type (e.g., "CD4Naive")
chr <- args[3]          # Chromosome (e.g., "chr1")
intron <- args[4]       # Intron/junction cluster ID to test

# Determine trait type (case-control vs. quantitative)
# Affects COLOC input parameters and priors
num_1 <- which(cc_trait == gwas_name)
num_2 <- which(quant_trait == gwas_name)

if(length(num_1) > 0) trait_type <- "cc"       # Case-control
if(length(num_2) > 0) trait_type <- "quant"    # Quantitative

# Get cell type sample size index
num_3 <- which(cell_type == celltype)

# Load required libraries
library("coloc")         # Colocalization analysis (Giambartolomei et al. 2014)
library(dplyr)           # Data manipulation
library(data.table)      # Fast file reading

#Load and Prepare Input Data
##=============================================================================
## LOAD GWAS SUMMARY STATISTICS
##=============================================================================
# Read GWAS summary statistics for specific locus
# Data organized by: celltype/chromosome/gwas_trait/intron_cluster
# This structure enables efficient parallel processing

GWAS_data <- paste0(celltype, "_", chr, "_nominals_1/", gwas_name, "/", intron)
print("start read in")

# Read GWAS data
# Expected columns: variant_id, beta, se, pval_nominal
GWAS_pre <- fread(GWAS_data, sep = "\t", header = TRUE)
GWAS_pre <- as.data.frame(GWAS_pre)
GWAS_pre <- GWAS_pre[, 2:5]  # Keep only essential columns
colnames(GWAS_pre) <- c("variant_id", "beta", "se", "pval_nominal")

##=============================================================================
## LOAD MINOR ALLELE FREQUENCY DATA
##=============================================================================
# MAF required for COLOC approximate Bayes factor calculation
# Helps estimate variance of effect sizes
# Should match reference panel used for GWAS

MAF <- read.table(paste0(chr, ".txt"), header = TRUE)
colnames(MAF) <- c("variant_id", "maf")

print("finished prepare")

##=============================================================================
## INITIALIZE RESULTS DATA FRAME
##=============================================================================
# Prepare output structure for colocalization results
# 16 columns capture key information about shared signals

test_rs <- as.data.frame(matrix(NA, 1, 16))
colnames(test_rs) <- c(
    "sqtl_variant_id_most",    # Top sQTL variant (before merging)
    "intron_id",               # Junction cluster ID
    "gene",                    # Gene symbol
    "nsnps",                   # Number of SNPs tested in region
    "PP.H0.abf",              # Posterior prob: no association with either trait
    "PP.H1.abf",              # Posterior prob: association with trait 1 only (GWAS)
    "PP.H2.abf",              # Posterior prob: association with trait 2 only (sQTL)
    "PP.H3.abf",              # Posterior prob: both traits, different causal variants
    "PP.H4.abf",              # Posterior prob: both traits, SHARED causal variant
    "sqtl_variant_P_most",    # P-value of top sQTL (before merging)
    "sqtl_variant_id_merge",  # Top sQTL variant (after merging with GWAS)
    "sqtl_variant_P_merge",   # P-value of top sQTL (after merging)
    "GWAS_variant_id_most",   # Top GWAS variant (before merging)
    "GWAS_variant_P_most",    # P-value of top GWAS variant (before merging)
    "GWAS_variant_id_merge",  # Top GWAS variant (after merging with sQTL)
    "GWAS_variant_P_merge"    # P-value of top GWAS variant (after merging)
)

#Extract Top Variants and Merge Datasets
##=============================================================================
## LOAD sQTL SUMMARY STATISTICS
##=============================================================================
# Read sQTL nominal pass results for this intron cluster
# Contains all variant-junction associations within cis-window (1 Mb)

snp_1mb <- fread(paste0(celltype, "_", chr, "_nominals_1/", gwas_name, "/", 
                        intron), sep = " ", header = FALSE)
colnames(snp_1mb) <- c("gene", "intron", "variant_id", "pval_nominal", "beta")
snp_1mb <- as.data.frame(snp_1mb)

# Sort by P-value to identify top sQTL variant
snp_1mb <- snp_1mb[order(as.numeric(snp_1mb$pval_nominal)), ]

##=============================================================================
## RECORD TOP VARIANTS BEFORE MERGING
##=============================================================================
# Document lead variants from each dataset independently
# Important for tracking potential differences after filtering

# Top sQTL variant (most significant in sQTL analysis)
test_rs[1, 1] <- snp_1mb[1, 3]      # variant_id
test_rs[1, 2] <- snp_1mb[1, 2]      # intron_id
test_rs[1, 3] <- snp_1mb[1, 1]      # gene
test_rs[1, 10] <- snp_1mb[1, 4]     # P-value

# Top GWAS variant (most significant in GWAS)
GWAS_pre_order <- GWAS_pre[order(as.numeric(GWAS_pre$pval_nominal)), ]
test_rs[1, 13] <- GWAS_pre_order[1, 1]  # variant_id
test_rs[1, 14] <- GWAS_pre_order[1, 4]  # P-value

##=============================================================================
## MERGE DATASETS FOR COLOCALIZATION
##=============================================================================
# Create unified dataset containing both sQTL and GWAS statistics
# Only variants present in BOTH datasets are retained
# This is critical for proper colocalization testing

# Step 1: Merge sQTL with MAF data
input <- merge(snp_1mb, MAF, by = "variant_id")

# Step 2: Merge with GWAS data
# suffixes distinguish columns from each dataset
input <- merge(input, GWAS_pre, by = "variant_id", 
               suffixes = c("_sqtl", "_gwas"))

# Sort by GWAS P-value for reporting
input <- input[order(as.numeric(input$pval_nominal_gwas)), ]

##=============================================================================
## RECORD TOP VARIANTS AFTER MERGING
##=============================================================================

# Top sQTL variant (in merged dataset)
input_ordersqtl <- input[order(as.numeric(input$pval_nominal_sqtl)), ]
test_rs[1, 11] <- input_ordersqtl[1, 1]  # variant_id
test_rs[1, 12] <- input_ordersqtl[1, 4]  # P-value

# Top GWAS variant (in merged dataset)
input_ordergwas <- input[order(as.numeric(input$pval_nominal_gwas)), ]
test_rs[1, 15] <- input_ordergwas[1, 1]  # variant_id
test_rs[1, 16] <- input_ordergwas[1, 9]  # P-value

#Run COLOC Analysis
##=============================================================================
## CALCULATE VARIANCE OF BETA (REQUIRED FOR COLOC)
##=============================================================================
# Variance of effect size = (standard error)^2
# Used in approximate Bayes factor calculation
input$varbeta <- (input$se)^2

##=============================================================================
## PERFORM COLOCALIZATION ANALYSIS
##=============================================================================
# COLOC uses approximate Bayes factors to assess evidence for shared causal variants
# Tests 5 hypotheses (H0-H4) and calculates posterior probabilities
#
# Hypotheses:
# H0: No association with either trait
# H1: Association with GWAS only (no sQTL)
# H2: Association with sQTL only (no GWAS)
# H3: Both traits associated, but different causal variants (linkage)
# H4: Both traits associated, SHARED causal variant (colocalization)
#

# Minimum variant threshold: At least 10 SNPs required for reliable analysis
if(nrow(input) > 10) {
    
    ##=========================================================================
    ## CASE-CONTROL GWAS
    ##=========================================================================
    if(trait_type == "cc") {
        result <- coloc.abf(
            # Dataset 1: GWAS (case-control)
            dataset1 = list(
                pvalues = input$pval_nominal_gwas,  # P-values
                type = "cc",                         # Case-control trait
                beta = input$beta_gwas,              # Effect sizes
                varbeta = input$varbeta,             # Variance of beta
                s = cc_trait_case[num_1] / (cc_trait_case[num_1] + 
                                            cc_trait_control[num_1]),
                # s = case proportion (affects priors)
                N = cc_trait_case[num_1] + cc_trait_control[num_1],
                # Total sample size
                snp = input$variant_id               # Variant IDs
            ),
            
            # Dataset 2: sQTL (always quantitative - normalized PSI values)
            dataset2 = list(
                pvalues = input$pval_nominal_sqtl,  # P-values
                type = "quant",                      # Quantitative trait
                N = celltype_sample[num_3],          # Sample size for cell type
                snp = input$variant_id               # Variant IDs
            ),
            
            MAF = input$maf  # Minor allele frequencies (improves accuracy)
        )
        
        # Extract posterior probabilities for 5 hypotheses + nsnps
        test_rs[1, 4:9] <- t(as.data.frame(result$summary))[1, 1:6]
    }
    
    ##=========================================================================
    ## QUANTITATIVE GWAS
    ##=========================================================================
    if(trait_type == "quant") {
        result <- coloc.abf(
            # Dataset 1: GWAS (quantitative)
            dataset1 = list(
                pvalues = input$pval_nominal_gwas,  # P-values
                type = "quant",                      # Quantitative trait
                beta = input$beta_gwas,              # Effect sizes
                varbeta = input$varbeta,             # Variance of beta
                N = quant_trait_num[num_2],          # Sample size
                snp = input$variant_id               # Variant IDs
            ),
            
            # Dataset 2: sQTL (quantitative)
            dataset2 = list(
                pvalues = input$pval_nominal_sqtl,  # P-values
                type = "quant",                      # Quantitative trait
                N = celltype_sample[num_3],          # Cell type sample size
                snp = input$variant_id               # Variant IDs
            ),
            
            MAF = input$maf  # Minor allele frequencies
        )
        
        # Extract posterior probabilities
        test_rs[1, 4:9] <- t(as.data.frame(result$summary))[1, 1:6]
    }
}


#Write Results
##=============================================================================
## WRITE OUTPUT
##=============================================================================
# Save colocalization results for this locus
# Output organized by cell type, chromosome, and GWAS trait

if(nrow(test_rs) > 0) {
    # Create output directory structure
    output_dir <- paste0(celltype, "_", chr, "_nominals_1/", gwas_name)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Write results
    # File named by intron cluster ID for easy lookup
    write.table(test_rs, 
                file = paste0(output_dir, "/", intron), 
                sep = "\t", 
                quote = FALSE, 
                row.names = FALSE, 
                col.names = FALSE)
}