================================================================================
Doublet Detection and Removal - Input/Output Description
================================================================================

================================================================================
INPUT FILE
================================================================================
{sample}_raw.rds

Format   : R binary file (RDS), Seurat object
Content  : Preprocessed single-cell RNA-seq data for one sample

Origin   : Generated from 10x Chromium sequencing data through the
following steps:
1. Loaded into R with Read10X, converted to Seurat object
2. Randomly subsampled to 1,000 cells
[Note: this step is only for testing; not required in production]
3. Preprocessed: normalize -> HVG (3,000) -> scale -> PCA (30 PCs) -> UMAP -> clustering
[Note: detailed preprocessing commands are recorded in the {sample}_raw.rds object and can be retrieved via seurat_obj@commands]

Contains :

Raw UMI count matrix
Normalized expression data
PCA and UMAP embeddings

Scale: 1,000 cells (subsampled for testing), all genes

================================================================================
OUTPUT FILES
================================================================================
{sample}_doublet_removed.rds

Format   : R binary file (RDS), Seurat object
Content  : Input Seurat object with two additional metadata columns:

DF.classfication1   Initial doublet classification (unadjusted nExp)
DF.classfication2   Final doublet classification (homotypic-adjusted nExp)

Values: Singlet / Doublet

Usage: Filter doublets for downstream analysis:
subset(seurat_obj, subset = DF.classfication2 == "Singlet")

{sample}_meta.tsv

Format   : Tab-delimited text, rows = cells, columns = metadata fields
Content  : Per-cell metadata table exported from the Seurat object
Used for QC reporting and doublet rate statistics