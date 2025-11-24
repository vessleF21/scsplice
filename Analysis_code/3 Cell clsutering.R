
##Harmony Integration and Cell Type Annotation PipelineCommand-line Interface Setup###

# Official Seurat repository link: https://satijalab.org/seurat/
# Official harmony repository link: https://github.com/immunogenomics/harmony

# Parse command-line arguments for batch integration and clustering

library(optparse)
option_list = list(
  make_option("--filename", action="store", default=NA, type='character', 
              help="Path to input Seurat object (.qs file) [required]"),
  make_option("--batchss", action="store", default=NA, type='character', 
              help="Batch variable name for integration (e.g., 'POOL') [required]"),
  make_option("--dimss", action="store", default=0, type='integer', 
              help="Number of dimensions for downstream analysis [required]"),
  make_option("--cellTY", action="store", default=NA, type='character', 
              help="Cell type label for output file naming [required]"),
  make_option("--filter", action="store_true", default=FALSE, 
              help="Apply additional filtering if TRUE"),
  make_option("--neighbors_within_batch", action="store_true", default=FALSE, 
              help="Compute neighbors within batch if TRUE")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load required libraries
library(Seurat)
library(patchwork)
library(harmony)        # For batch effect correction
library(tidyverse)
library(clustree)
library(cowplot)
library(qs)            # Fast serialization for large Seurat objects
library(bbknnR)
library(dplyr)
library(stringr)
library(ggplot2)
options(future.globals.maxSize = 1000 * 1024^3)  # Increase memory limit for large datasets

# Load aggregated Seurat object from all pools
sce = qread(file = opt$filename)

# Remove previous clustering results to avoid confusion
# Cleans up resolution-specific clustering columns from prior analyses
index <- match(paste0("RNA_snn_res.", seq(0.1, 3, by=0.1)), 
               colnames(sce@meta.data))
index <- index[!is.na(index)]
sce@meta.data <- sce@meta.data[,-index]

# Data Normalization and Feature Selection
# Log-normalization as specified in methodology:
# "We performed log-normalization using the NormalizeData function"
# Scale factor of 10,000 is standard for 10x Chromium data
sce <- NormalizeData(object = sce, 
                     normalization.method = "LogNormalize", 
                     scale.factor = 1e4)

# Identify highly variable features for dimensionality reduction
# Focuses analysis on genes with greatest variation across cells
sce <- FindVariableFeatures(sce)

# Scale data using linear regression as mentioned in methodology:
# "Linear regression using ScaleData functions"
sce <- ScaleData(object = sce, features = rownames(sce))

# PCA for initial dimensionality reduction
# Uses variable features to capture major sources of variation
# 50 PCs computed to allow flexible downstream analysis
sce <- RunPCA(object = sce, 
              features = VariableFeatures(object = sce), 
              npcs = 50)  # Apply Harmony for batch effect correction as specified in methodology:
# "We applied the RunHarmony function in harmony (v1.2.0) to reduce 
#  dimensionality and mitigate batch effects"
# 
# Corrects for technical variation between pools while preserving biological signal
# Each pool was processed independently, creating batch effects that need correction
sce <- sce %>% 
  RunHarmony(group.by.vars = "POOL")  
# POOL = batch/pool identifier from opt$batchss Neighbor Finding and UMAP Dimensionality Reductionr
# Find k-nearest neighbors in Harmony-corrected space
# Use user-specified number of dimensions (opt$dimss)
# This addresses: "high-dimensional variables are common in scRNA-seq data"
sce <- sce %>% 
  FindNeighbors(reduction = "harmony", dims = 1:opt$dimss) %>% 
  RunUMAP(reduction = "harmony", dims = 1:opt$dimss, verbose = F)
# Generate QC plots to assess batch integration quality
# p1: Check if batch effects are removed (pools should be well-mixed)
# p2: Visualize Level 2 cell type annotations from OneK1K reference
# p3: Visualize detailed cell type annotations from OneK1K reference
p1 <- DimPlot(sce, reduction = "umap", group.by = "POOL", raster = F)
p2 <- DimPlot(sce, reduction = "umap", group.by = "l2_onek1k", 
              label = TRUE, raster = F)
p3 <- DimPlot(sce, reduction = "umap", group.by = "cell_type_onek1k", 
              label = TRUE, raster = F)

# Combine plots for comprehensive visualization
P.total = wrap_plots(p1, p2, p3)
ggsave(plot = P.total, 
       filename = paste0(opt$cellTY, "_", opt$dimss, "_", opt$batchss, 
                        "_Step6.Umap_plot.png"),
       width = 48, height = 6)

# Save intermediate object before clustering
qsave(sce, file = paste0(opt$cellTY, "_", opt$dimss, "_", opt$batchss, 
                         "_before_sce.harmony.qs"))  # Perform clustering at multiple resolutions as mentioned in methodology:
# "Cell clusters were identified using the FindClusters function in Seurat"

# Multiple resolutions allow identification of both major cell types and subtypes
for (res in c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)) {
  print(res)
  sce <- FindClusters(sce, 
                      resolution = res, 
                      algorithm = 1)  # Louvain algorithm (default)
}

# Save object with all clustering resolutions
qsave(sce, file = paste0(opt$cellTY, "_", opt$dimss, "_", opt$batchss, 
                         "_batch_cov_sce_qc_after_findcluster.qs"))
                         # Generate comprehensive visualization of all clustering resolutions
# Allows visual assessment of optimal resolution for cell type annotation
# Methodology states clusters are annotated based on marker gene expression:
#   - B cells: CD19, MS4A1
#   - Plasma cells (PB): JCHAIN, MZB1
#   - T cells: CD3D, CD3E, IL7R
#   - NK cells: NKG7, KLRF1
#   - Monocytes: CD14, FCGR3A
#   - DC: CD1C, CLEC10A
#   - Progenitors: SOX4
#   - Platelets: PPBP, PF4

cluster_umap <- wrap_plots(ncol = 4,
  DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.3", 
          label = T, raster=FALSE) & NoAxes(),
  DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.4", 
          label = T, raster=FALSE) & NoAxes(),
  DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.5", 
          label = T, raster=FALSE) & NoAxes(), 
  DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.6", 
          label = T, raster=FALSE) & NoAxes(),
  DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.7", 
          label = T, raster=FALSE) & NoAxes(),
  DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.8", 
          label = T, raster=FALSE) & NoAxes(),
  DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.9", 
          label = T, raster=FALSE) & NoAxes(), 
  DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.1", 
          label = T, raster=FALSE) & NoAxes(),
  DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.1.1", 
          label = T, raster=FALSE) & NoAxes(), 
  DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.1.2", 
          label = T, raster=FALSE) & NoAxes()
)

ggsave(cluster_umap,
       filename = paste0(opt$cellTY, "_", opt$dimss, "_", opt$batchss, 
                        "_batch_cov", "_BBKNN_Umap_plot_2.png"),
       width = 30, height = 10)