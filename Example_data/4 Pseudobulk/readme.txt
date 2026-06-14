================================================================================
Pseudobulk BAM Generation - Input/Output Description
================================================================================

================================================================================
OVERVIEW
================================================================================
Splits a BAM file into per-cell-type BAM files based on cell barcode
annotations. Filters reads by barcode whitelist and removes multimappers.

================================================================================
INPUT FILES
================================================================================

--------------------------------------------------------------------------------
1. {sample}_100pct.bam
--------------------------------------------------------------------------------
Format   : BAM (sorted, indexed)
Content  : Subsampled (1%) aligned reads from STARsolo output
           Deduplicated with Picard MarkDuplicates

Required tags per read:
  CB    Cell barcode tag (assigned by STARsolo)
  NH    Number of alignments (for multimapper filtering)
  vW    WASP filter tag (optional; if absent, NH==1 filter is used)

--------------------------------------------------------------------------------
2. {sample}@{celltype}.txt
--------------------------------------------------------------------------------
Format   : Tab-delimited text, no header
Content  : Cell barcode to cell type mapping for one cell type

Example  :
  AACCGCGTCTTGGGTA    B1@682_683@CD4Naive
  AACGTTGTCGAATGCT    B1@682_683@CD4Naive
  AACTCAGCACGACGAA    B1@682_683@CD4Naive

Columns  :
  Col 1   Cell barcode (must match CB tag in BAM)
  Col 2   Cell type label

Origin   : Generated from Seurat cell metadata after cell type annotation,
           split per sample per cell type

================================================================================
OUTPUT FILES
================================================================================

--------------------------------------------------------------------------------
{out_prefix}.{celltype}.bam
--------------------------------------------------------------------------------
Format   : BAM (unsorted, requires indexing after creation)
Content  : Reads from cells belonging to the specified cell type only

Filtering applied:
  - Reads with CB tag not in barcode whitelist are discarded
  - If vW tag present: only reads with vW==1 (WASP pass) are retained
  - If vW tag absent : only uniquely mapped reads (NH==1) are retained
  - Multimappers (NH>1) without vW tag are discarded

--------------------------------------------------------------------------------
{out_prefix}.{celltype}.bam.bai
--------------------------------------------------------------------------------
Format   : BAM index
Content  : Index file for the output BAM, required for downstream tools
           (e.g. regtools junctions extract, IGV visualization)


