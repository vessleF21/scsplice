



# =============================================================================
# Part A. BAM File Splitting by Sample
# Filter and split BAM file by sample-specific barcodes
# Extracts reads belonging to singlet cells assigned to single sample-onek1k_sample1
# onek1k_sample1.txt: barcode list extracted from Demuxlet output
# Uses 20 threads for parallel processing
# Official sinto repository link: https://github.com/timoast/sinto
# =============================================================================

sinto filterbarcodes -p 20 \
                     -b possorted_genome_bam.bam \
                     -c onek1k_sample1.txt \
                     --outdir ./split_bam/onek1k_sample1


# =============================================================================
# Part B. VCF File Processing for Individual Samples
# Extract sample-specific variants with at least 1 alternate allele
# Removes monomorphic sites for the specific sample
# Official bcftools repository link: https://github.com/samtools/bcftools
# =============================================================================

# J1. Extract sample-specific variants with at least 1 alternate allele
bcftools view -c1 -s onek1k_sample1 onek1k_sample1_vcf

# J2. Sort VCF file and compress with bgzip
bcftools sort -Oz -o onek1k_sample1.vcf.gz

# J3. Create index for compressed VCF
bcftools index onek1k_sample1.vcf.gz


# =============================================================================
# Part C. STAR Alignment with WASP Filtering
# Two-pass STAR alignment with single-cell and allele-specific mapping features
# Key parameters:
# - Uses GRCh38 reference genome
# - Forward strand specificity for 10x Chromium chemistry
# - CB_UMI_Simple mode for 10x data processing
# - 1MM_Directional_UMItools for UMI deduplication
# - WASP filtering to remove reference allele mapping bias
# - Outputs SAM attributes required for downstream QC and quantification
# Official STAR repository link: https://github.com/alexdobin/STAR
# =============================================================================

STAR \
    --runThreadN 24 \
    --genomeDir ./star/hg38_gencode_v32 \
    --twopassMode Basic \
    --soloStrand Forward \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist ./ref/737K-august-2016.txt \
    --soloCBmatchWLtype 1MM \
    --soloUMIdedup 1MM_Directional_UMItools \
    --readFilesIn $readFilesIn \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM XS vG vA vW \
    --waspOutputMode SAMtag \
    --varVCFfile <(zcat onek1k_sample1.vcf.gz) \
    --sjdbGTFfile ./hg38/gencode.annotation.gtf \
    --soloOutFileNames onek1k_sample1 \
    --outFileNamePrefix onek1k_sample1


# =============================================================================
# Part D. PCR Duplicate Removal
# BARCODE_TAG=CB: Uses cell barcode to identify duplicates within each cell
# MOLECULAR_IDENTIFIER_TAG=UB: Uses UMI-corrected barcode for molecular duplicate detection
# REMOVE_DUPLICATES=true: Physically removes duplicates from output BAM
# Official STAR repository link: https://broadinstitute.github.io/picard/
# =============================================================================

java -jar picard.jar MarkDuplicates \
    REMOVE_DUPLICATES=true \
    BARCODE_TAG=CB \
    MOLECULAR_IDENTIFIER_TAG=UB \
    I=onek1k_sample1.sorted.bam \
    O=onek1k_sample1.sorted_rmdup.bam \
    M=marked_dup_metrics.onek1k_sample1.txt


##############################################################################
## Alternative Splicing (AS) Quantification Pipeline
## Official regtools repository link: https://github.com/griffithlab/regtools
## Official leafcutter repository link: https://github.com/davidaknowles/leafcutter
##############################################################################

# =============================================================================
# Part E. Splice Junction Extraction with RegTools
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
# =============================================================================

regtools junctions extract \
    -a 5 \           # Anchor length: minimum overlap on each exon
    -m 30 \          # Min intron size: filters short junctions
    -M 500000 \      # Max intron size: filters ultra-long junctions
    -s XS \          # Strand tag from STAR alignment
    bam \            # Input: pseudobulk BAM (cell type-specific, WASP-filtered)
    -o bam.junc      # Output: junction file with coordinates and counts


# =============================================================================
# Part F. Junction Clustering with LeafCutter
# Cluster junctions that share splice sites using LeafCutter
# Implements: "junctions that shared splice sites were clustered using 
# leafcutter_cluster_regtools.py script"
#
# This groups alternative splice junctions into intron clusters
#
# Parameters:
#   -j: Input file listing all junction files (one per sample)
#   -m 20: Minimum reads supporting a cluster across all samples
#   -o: Output prefix for cluster files
#   -l 500000: Maximum intron length
# =============================================================================

python leafcutter_cluster_regtools.py \
    -j test_juncfiles.txt_celltype \  # List of junction files per sample
    -m 20 \                            # Min reads per cluster
    -o oneK1K_celltype \               # Output prefix
    -l 500000                          # Max intron length filter

# Output files:
#   - oneK1K_celltype_perind.counts.gz: Per-individual junction counts


# =============================================================================
# Part G. Normalize Junction Usage
# Calculate normalized splice junction usage (Percent Spliced In - PSI)
# Implements: "normalized splice junction usage across samples was calculated 
# using prepare_phenotype_table.py script"
#
# This script:
# 1. Converts raw counts to proportions within each intron cluster
# 2. Normalizes for sequencing depth
# 3. Prepares data for QTL mapping
# =============================================================================

python prepare_phenotype_table.py \
    oneK1K_celltype_perind.counts.gz


# =============================================================================
# Part H. Junction Quality Filtering
# Filter junctions based on detection rate and variability
# Implements methodology criteria:
# 1. "detected in < 40% of samples" -> EXCLUDE
# 2. "minimal variability (standard deviation < 0.005)" -> EXCLUDE
#
# These filters ensure:
#   - Robust detection across cohort
#   - Sufficient variability for association testing
# =============================================================================

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
# Format: One intron ID per line
write.table(
    introns, 
    paste0('filtered_celltype.txt'), 
    col.names=F, 
    row.names=F, 
    quote=F
)


# =============================================================================
# Part I. Monocyte Trajectory AS Analysis
# Dynamic AS during cell differentiation: Monocyte trajectory analysis
# PHATE dimensionality reduction -> Slingshot pseudotime -> Q4 quantile binning
# -> Long-format PSI matrix -> One-way ANOVA per intron junction
#
# Input : mono.data.qs, perind_counts.tsv
# Output: metadata.tsv, tmp_long.txt, anova_results.csv
# =============================================================================
library(Seurat)
library(qs)
library(dplyr)
library(reticulate)
library(phateR)
library(slingshot)
library(data.table)

options(future.globals.maxSize = 1000 * 1024^3)

# I1. PHATE -> Slingshot pseudotime -> Q4 binning
# PHATE reduces high-dimensional gene expression to 2D embedding for trajectory
# inference. Slingshot fits a principal curve through CD14 monocytes as the
# start cluster to assign pseudotime. Cells are then binned into Q1-Q4 quantiles
# for downstream differential splicing analysis.

# Load monocyte Seurat object from compressed file
obj <- qread("mono.data.qs")

# Set RNA as the active assay for downstream analysis
DefaultAssay(obj) <- "RNA"

# Run PHATE dimensionality reduction on scaled gene expression data
# npca=10: first reduce to 10 PCs before PHATE embedding
# Store result as a new dimensional reduction slot "phate" in Seurat object
obj[["phate"]] <- CreateDimReducObject(
  phate(t(GetAssayData(obj, slot="scale.data")), npca=10)$embedding,
  key="PHATE_", assay="RNA")

# Fit trajectory using Slingshot on PHATE embedding
# clusterLabels: use cell type annotations to guide trajectory
# start.clus="CD14": set CD14 monocytes as the root of the trajectory
traj <- slingshot(Embeddings(obj,"phate"), clusterLabels=obj$CT_res, start.clus="CD14")

# Extract pseudotime values (arc length along principal curve) and add to Seurat metadata
obj  <- AddMetaData(obj, setNames(as.data.frame(slingCurves(traj)[[1]]$lambda), "pseudotime"))

# Copy barcode column from pool_sample_Barcode for downstream joining
obj@meta.data$barcode <- obj@meta.data$pool_sample_Barcode

# Join pseudotime values (Lineage1) to all cells using barcode as key
# full_join ensures all cells are retained even if pseudotime is NA
pt   <- dplyr::full_join(
          dplyr::select(tibble::as_tibble(obj[[]]), barcode),
          dplyr::select(tibble::rownames_to_column(
            as.data.frame(slingPseudotime(traj, na=FALSE)), "barcode"),
            barcode, pt=Lineage1), by="barcode")

# Bin pseudotime into 4 quantile groups (Q1-Q4) using equal-frequency binning
pt$Q4 <- `levels<-`(Hmisc::cut2(pt$pt, g=4), paste0("Q",1:4))

# Add Q4 bin labels back to Seurat object metadata
obj   <- AddMetaData(obj, tibble::column_to_rownames(
           dplyr::select(pt, barcode, Q4), "barcode"))

# Export full metadata table including pseudotime and Q4 bins
write.table(as.data.frame(obj@meta.data), "metadata.tsv",
            sep="\t", row.names=TRUE, quote=FALSE)


# I2. Preprocess perind_counts.tsv into long format (awk)
# Parses wide-format junction counts (numerator/denominator) into long format.
# Skips zero-coverage cells (0/0); computes PSI = numerator/denominator.
# Column header format: individual@batch@bin; bin (Q1-Q4) is extracted as group label.
# Output: tmp_long.txt (intron | PSI | bin)
system("awk 'NR==1{for(i=2;i<=NF;i++){split($i,h,\"@\"); hdr[i]=h[3]}; next}
             {for(i=2;i<=NF;i++){if($i==\"0/0\") next;
              split($i,v,\"/\"); print $1\"\\t\"v[1]/v[2]\"\\t\"hdr[i]}}' \
             perind_counts.tsv > tmp_long.txt")

# I3. One-way ANOVA per intron junction
# Tests whether PSI differs significantly across pseudotime quantiles (Q1-Q4).
# Model: PSI ~ bin, run independently per intron. Introns where the model
# fails (e.g. only one group present) are dropped via na.omit.
# P values adjusted using the Holm method.
# Output: anova_results.csv (intron | F | p | p_adj)
tbl <- fread("tmp_long.txt", header=FALSE, col.names=c("intron","psi","bin"))
out <- na.omit(do.call(rbind, lapply(unique(tbl$intron), function(jid) {
  res <- summary(aov(psi ~ bin, data=tbl[intron==jid]))[[1]]
  data.frame(intron=jid, F=res["bin","F value"], p=res["bin","Pr(>F)"])
})))
out$p_adj <- p.adjust(out$p)
write.csv(out, "anova_results.csv", row.names=FALSE)


# =============================================================================
# Part J. Demographic-Associated Differential Splicing Analysis
# Identify sex- and age-associated differential splicing events
# Implements: "We identified demographic-associated differential splicing events 
# using the leafcutter_ds.R script"
#
# Analysis strategy:
# - Sex: Compare male vs. female
# - Age: Pairwise comparisons between three groups:
#     * Younger (<=45 years)
#     * Middle-aged (45-60 years)
#     * Older (>=60 years)
# - Statistical significance: FDR < 0.05 (Benjamini-Hochberg correction)
# =============================================================================

# Define paths
Rscript_path="R/bin/Rscript"
script_path="leafcutter_ds.R"       # LeafCutter differential splicing script
output_dir="./result/"
input_dir="./junc/"
ct_dir="./ct/" 

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