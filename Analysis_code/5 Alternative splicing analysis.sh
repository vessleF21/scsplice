##Alternative Splicing (AS) Quantification Pipeline##

# Official regtools repository link: https://github.com/griffithlab/regtools
# Official leafcutter repository link: https://github.com/davidaknowles/leafcutter


#Step 1: Splice Junction Extraction with RegTools
# Extract splice junctions from pseudobulk BAM files using RegTools
# This implements: "exon-exon junctions were identified and extracted based on 
# CIGAR strings from pseudobulk BAM files"
#
# Parameters aligned with methodology's filtering criteria:
#   -a 5: Minimum anchor length (5 bp on each side of junction)
#   -m 30: Minimum intron length (30 bp) - filters out very short junctions
#   -M 500000: Maximum intron length (500 kb) - filters out extremely long junctions
#   -s XS: Use XS tag for strand information (from STAR alignment)
#
# Quality filters applied (from methodology):
# 1. >= 5 supporting reads (enforced by RegTools default or downstream filtering)
# 2. >= 30 bp junction length (-m 30)
# 3. <= 500 kb junction length (-M 500000)

regtools junctions extract \
    -a 5 \           # Anchor length: minimum overlap on each exon
    -m 30 \          # Min intron size: filters short junctions
    -M 500000 \      # Max intron size: filters ultra-long junctions
    -s XS \          # Strand tag from STAR alignment
    bam \            # Input: pseudobulk BAM (cell type-specific, WASP-filtered)
    -o bam.junc      # Output: junction file with coordinates and counts


#Step 2: Junction Clustering with LeafCutter
# Cluster junctions that share splice sites using LeafCutter
# Implements: "junctions that shared splice sites were clustered using 
# leafcutter_cluster_regtools.py script"
#
# This groups alternative splice junctions into intron clusters
# Example: Three isoforms sharing donor/acceptor sites form one cluster
#
# Parameters:
#   -j: Input file listing all junction files (one per sample)
#   -m 20: Minimum reads supporting a cluster across all samples
#   -o: Output prefix for cluster files
#   -l 500000: Maximum intron length (consistent with RegTools filtering)

python leafcutter_cluster_regtools.py \
    -j test_juncfiles.txt_celltype \  # List of junction files per sample
    -m 20 \                            # Min reads per cluster
    -o oneK1K_celltype \               # Output prefix
    -l 500000                          # Max intron length filter

# Output files:
#   - oneK1K_celltype_perind.counts.gz: Per-individual junction counts
#   - oneK1K_celltype_pooled: Pooled junction information
#Step 3: Normalize Junction Usage
# Calculate normalized splice junction usage (Percent Spliced In - PSI)
# Implements: "normalized splice junction usage across samples was calculated 
# using modified prepare_phenotype_table.py script"
#
# This script:
# 1. Converts raw counts to proportions within each intron cluster
# 2. Normalizes for sequencing depth
# 3. Prepares data for QTL mapping (phenotype matrix)

python prepare_phenotype_table.py \
    oneK1K_celltype_perind.counts.gz

# Output:
#   - oneK1K_celltype_perind_numers.counts.gz: Numerator counts for PSI
#   - oneK1K_celltype_perind_denoms.counts.gz: Denominator counts for PSI

#Step 4: Junction Quality Filtering
##=============================================================================
## Intron Quality Filtering
##=============================================================================
# Filter junctions based on detection rate and variability
# Implements methodology criteria:
# 1. "detected in < 40% of samples" -> EXCLUDE
# 2. "minimal variability (standard deviation < 0.005)" -> EXCLUDE
#
# These filters ensure:
#   - Robust detection across cohort
#   - Sufficient variability for association testing

# Parse command-line argument: cell type identifier
args <- commandArgs(trailingOnly = TRUE)
ct <- args[1]

# Load normalized junction counts (numerators for PSI calculation)
# Each row = intron/junction
# Each column = sample
dat <- read.delim(
    paste0('oneK1K_celltype_perind_numers.counts.gz'), 
    sep=' ', 
    check.names=F
)

# Initialize vector to store passing introns
introns <- c()
colNum <- ncol(dat)

# Filter based on detection rate
# Criterion: Junction must be detected (non-zero) in >= 40% of samples
# This implements: "detected in < 40% of samples were filtered out"
for (i in 1:nrow(dat)){
    # Calculate proportion of samples with zero counts
    zero_proportion <- length(which(dat[i,] == 0)) / colNum
    
    # Retain if detected in >= 40% of samples (i.e., zeros in <= 60%)
    if (zero_proportion <= 0.6){
        introns <- c(introns, rownames(dat)[i])
    }
}

# Save filtered intron list for downstream sQTL analysis
# Format: One intron ID per line (chr:start:end:cluster_id)
write.table(
    introns, 
    paste0('filtered_celltype.txt'), 
    col.names=F, 
    row.names=F, 
    quote=F
)

# Note: The standard deviation filter (SD < 0.005) is typically applied
# within the LeafCutter differential splicing script or during sQTL analysis

#Step 5: Demographic-Associated Differential Splicing Analysis
##=============================================================================
## Demographic-Biased Alternative Splicing (AS) Detection
##=============================================================================
# Identify sex- and age-associated differential splicing events
# Implements: "We identified demographic-associated differential splicing events 
# using the leafcutter_ds.R script"
#
# Analysis strategy:
# - Sex: Compare male vs. female
# - Age: Pairwise comparisons between three groups:
#     * Younger (≤45 years)
#     * Middle-aged (45-60 years)  
#     * Older (≥60 years)
# - Statistical significance: FDR < 0.05 (Benjamini-Hochberg correction)

# Define paths
Rscript_path="R/bin/Rscript"
script_path="leafcutter_ds.R"       # LeafCutter differential splicing script
output_dir="./result/"
input_dir="./junc/"
ct_dir="./ct/"                      # Covariate files (sex, age groups)

# Loop through each cell type for demographic analysis
# This ensures cell type-specific differential splicing detection
for celltype in "${celltypes[@]}"; do
    echo "Processing $celltype..."
    
    # Output prefix for results
    output="${output_dir}leafcutter_ds.sex.${celltype}"
    
    # Input: Normalized junction counts for this cell type
    # Format: rows=junctions, columns=samples
    junc_file="${input_dir}S20_${celltype}/oneK1K_${celltype}_perind_numers.counts.gz"
    
    # Covariate file: Sample metadata with demographic information
    # Format: sample_id, sex, age_group, other_covariates
    # For sex analysis: binary variable (male/female)
    # For age analysis: categorical variable (younger/middle/older)
    ct_file="${ct_dir}${celltype}.txt"
    
    # Run LeafCutter differential splicing analysis
    # leafcutter_ds.R performs:
    # 1. Dirichlet-multinomial GLM for junction usage
    # 2. Likelihood ratio test for covariate effect
    # 3. FDR correction (Benjamini-Hochberg)
    # 4. Identification of significant differential splicing events (FDR < 0.05)
    $Rscript_path $script_path \
        -o $output \      # Output prefix
        $junc_file \      # Junction counts
        $ct_file          # Covariate file (sex/age)
    
    # Output files:
    #   - *.effect_sizes.txt: Effect sizes for each junction
    #   - *.cluster_significance.txt: Cluster-level significance
    #   - *.leafcutter_ds_cluster_summary.txt: Summary statistics
done
