================================================================================
Batch Correction, Clustering and Cell Type Annotation - Input/Output Description
================================================================================

================================================================================
INPUT FILE
================================================================================

--------------------------------------------------------------------------------
merged.qs
--------------------------------------------------------------------------------
Format   : Seurat object (.qs)
Content  : Merged single-cell RNA-seq data with batch (POOL) metadata

Origin   : Generated from Cell Ranger output through the following steps:
             1. Raw count matrix loaded with Read10X
             2. Split into 3 non-overlapping subsets to simulate batch structure
                [Note: for testing only; in production, load a pre-merged
                 multi-sample object with real POOL assignments]
             3. Merged into a single Seurat object
             4. POOL labels assigned: batch1 / batch2 / batch3

Contains :
  - Raw UMI count matrix
  - Cell metadata with POOL column (batch identity per cell)

Scale    : 3,000 cells (3 batches x 1,000 cells), all genes

================================================================================
OUTPUT FILES
================================================================================

--------------------------------------------------------------------------------
plot1.png
--------------------------------------------------------------------------------
Format   : PNG image
Content  : UMAP plot colored by POOL batch identity
           Used to assess whether Harmony successfully removed batch effects

--------------------------------------------------------------------------------
plot_2.png
--------------------------------------------------------------------------------
Format   : PNG image (30 x 10 inches)
Content  : 10-panel UMAP grid showing cluster assignments at each resolution
           (res 0.3 to 1.2, 4 panels per row)
           Used to visually select the optimal clustering resolution

--------------------------------------------------------------------------------
marker_dotplot.png
--------------------------------------------------------------------------------
Format   : PNG image (1200 x 600 px)
Content  : Dot plot of canonical marker genes across clusters at res 0.8
           Dot size = fraction of cells expressing the gene
           Dot color = average expression level
           Marker genes used:
             B cells      : CD19, MS4A1
             Plasma cells : JCHAIN, MZB1
             T cells      : CD3D, CD3E, IL7R
             NK cells     : NKG7, KLRF1
             Monocytes    : CD14, FCGR3A
             DC           : CD1C, CLEC10A
             Progenitors  : SOX4
             Platelets    : PPBP, PF4
