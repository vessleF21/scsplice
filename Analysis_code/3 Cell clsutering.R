# =============================================================================
# Batch correction, dimensionality reduction and clustering
# Uses Harmony to correct for batch effects across pools, then performs
# graph-based clustering at multiple resolutions for downstream cell type
# annotation.
#
# Input:
#   --filename              : path to input Seurat object (.qs)
#   --batchss               : batch variable name for UMAP labeling
#   --dimss                 : number of Harmony dimensions to use for
#                             neighbor graph and UMAP
#   --cellTY                : cell type label used in output filenames
#   --filter                : whether to apply additional filtering
#   --neighbors_within_batch: whether to restrict neighbors to same batch
#
# Output:
#   *_Step6.Umap_plot.png              : UMAP colored by pool, cell type labels
#   *_before_sce.harmony.qs            : Seurat object after Harmony, before clustering
#   *_batch_cov_sce_qc_after_findcluster.qs : Seurat object with cluster assignments
#   *_BBKNN_Umap_plot_2.png            : UMAP grid across resolutions 0.3-1.2
#
# Official links:
#   Harmony : https://github.com/immunogenomics/harmony
#   Seurat  : https://satijalab.org/seurat/
# =============================================================================
library(optparse)

option_list <- list(
  make_option("--filename",              type='character', help="Path to input .qs Seurat object"),
  make_option("--batchss",               type='character', help="Batch variable for UMAP labeling"),
  make_option("--dimss",                 type='integer',   help="Number of Harmony PCs to use"),
  make_option("--cellTY",                type='character', help="Cell type label for output filenames"),
  make_option("--filter",                action="store_true", default=FALSE, help="Apply additional filtering"),
  make_option("--neighbors_within_batch",action="store_true", default=FALSE, help="Restrict neighbors to same batch")
)
opt <- parse_args(OptionParser(option_list=option_list))

library(Seurat); library(patchwork); library(harmony)
library(tidyverse); library(clustree); library(cowplot)
library(qs); library(bbknnR); library(dplyr)
library(stringr); library(ggplot2)
options(future.globals.maxSize = 1000 * 1024^3)

# Load Seurat object
sce <- qread(opt$filename)

# Remove any existing clustering columns (RNA_snn_res.0.1 to 3.0)
# to avoid conflicts with downstream FindClusters calls
idx <- match(paste0("RNA_snn_res.", seq(0.1, 3, by=0.1)), colnames(sce@meta.data))
sce@meta.data <- sce@meta.data[, -idx[!is.na(idx)]]

# =============================================================================
# Step 1. Preprocessing: normalize, variable features, scale, PCA
# Standard Seurat pipeline before batch correction
# =============================================================================
sce <- NormalizeData(sce, normalization.method="LogNormalize", scale.factor=1e4)
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce, features=rownames(sce))
sce <- RunPCA(sce, features=VariableFeatures(sce), npcs=50)

# =============================================================================
# Step 2. Harmony batch correction
# Corrects PCA embeddings for batch effects across sequencing pools (POOL).
# Harmony iteratively clusters cells and adjusts embeddings to minimize
# differences between pools while preserving biological variation.
# =============================================================================
sce <- RunHarmony(sce, group.by.vars="POOL")

# =============================================================================
# Step 3. Neighbor graph and UMAP on Harmony embeddings
# Uses corrected Harmony dimensions (1:dimss) for both graph construction
# and UMAP visualization
# =============================================================================
sce <- FindNeighbors(sce, reduction="harmony", dims=1:opt$dimss)
sce <- RunUMAP(sce,      reduction="harmony", dims=1:opt$dimss, verbose=FALSE)

# =============================================================================
# Step 4. UMAP visualization before clustering
# Three panels: colored by pool, by OneK1K level-2 label, by OneK1K cell type
# =============================================================================
p1 <- DimPlot(sce, reduction="umap", group.by="POOL",              raster=FALSE)
p2 <- DimPlot(sce, reduction="umap", group.by="l2_onek1k",         label=TRUE, raster=FALSE)
p3 <- DimPlot(sce, reduction="umap", group.by="cell_type_onek1k",  label=TRUE, raster=FALSE)
ggsave(wrap_plots(p1,p2,p3),
       filename=paste0(opt$cellTY,"_",opt$dimss,"_",opt$batchss,"_Step6.Umap_plot.png"),
       width=48, height=6)

# Save Seurat object after Harmony correction, before clustering
qsave(sce, paste0(opt$cellTY,"_",opt$dimss,"_",opt$batchss,"_before_sce.harmony.qs"))

# =============================================================================
# Step 5. Graph-based clustering across multiple resolutions (0.3 - 1.2)
# Higher resolution -> more clusters. Multiple resolutions are tested to
# allow visual selection of the optimal granularity for cell type annotation.
# Algorithm 1 = Louvain
# =============================================================================
for (res in seq(0.3, 1.2, by=0.1))
  sce <- FindClusters(sce, resolution=res, algorithm=1)

qsave(sce, paste0(opt$cellTY,"_",opt$dimss,"_",opt$batchss,
                  "_batch_cov_sce_qc_after_findcluster.qs"))

# =============================================================================
# Step 6. UMAP grid across all resolutions for visual cluster inspection
# Each panel shows cluster assignments at one resolution; helps identify
# the resolution that best separates known cell types
# =============================================================================
res_list <- seq(0.3, 1.2, by=0.1)
plot_list <- lapply(res_list, function(res)
  DimPlot(sce, reduction="umap",
          group.by=paste0("RNA_snn_res.",res),
          label=TRUE, raster=FALSE) & NoAxes())

ggsave(wrap_plots(plot_list, ncol=4),
       filename=paste0(opt$cellTY,"_",opt$dimss,"_",opt$batchss,
                       "_batch_cov_BBKNN_Umap_plot_2.png"),
       width=30, height=10)