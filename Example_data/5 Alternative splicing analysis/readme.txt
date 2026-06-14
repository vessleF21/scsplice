================================================================================
LeafCutter Differential Splicing Analysis - Input/Output Description
================================================================================

================================================================================
OVERVIEW
================================================================================
Identifies differential splicing events between two groups (e.g. male vs female)
using LeafCutter. Pipeline proceeds from per-sample junction files through
intron clustering, phenotype normalization, and differential splicing analysis.

================================================================================
INPUT FILES
================================================================================

--------------------------------------------------------------------------------
1. {sample}@CD4Naive.bam.junc  (x6 samples)
   e.g. B1@682_683@CD4Naive.bam.junc
--------------------------------------------------------------------------------
Format   : BED-like tab-delimited text, no header
Content  : Splice junction reads extracted from cell type BAM files

Scale    : 6 samples (3 male, 3 female), cell type

--------------------------------------------------------------------------------
2. CD4Naive_juncfiles.txt
--------------------------------------------------------------------------------
Format   : Plain text, one file path per line
Content  : List of all 6 junction file paths

--------------------------------------------------------------------------------
3. ct_sex.txt
--------------------------------------------------------------------------------
Format   : Tab-delimited text, no header
Content  : Sample-to-group assignment for differential splicing analysis

Columns  :
  Col 1   Sample name (must match column names in counts.gz)
  Col 2   Group label (two groups only: male / female)

================================================================================
OUTPUT FILES
================================================================================

--------------------------------------------------------------------------------
oneK1K_celltype_perind_numers.counts.gz
--------------------------------------------------------------------------------
Format   : Gzip-compressed tab-delimited matrix
Content  : Raw intron usage counts per sample per intron cluster
           Rows = introns, Columns = samples

--------------------------------------------------------------------------------
oneK1K_celltype_perind.counts.gz
--------------------------------------------------------------------------------
Format   : Gzip-compressed tab-delimited matrix
Content  : Normalized PSI (percent spliced in) values per intron per sample

--------------------------------------------------------------------------------
oneK1K_celltype_perind.counts.gz.PCs
--------------------------------------------------------------------------------
Format   : Tab-delimited text
Content  : Top principal components of the normalized splicing matrix
           Used as covariates in QTL mapping to control for splicing
           heterogeneity

--------------------------------------------------------------------------------
CD4Naive_sex_cluster_significance.txt
--------------------------------------------------------------------------------
Format   : Tab-delimited text, one row per intron cluster
Content  : Differential splicing results per cluster

Columns  :
  cluster      Intron cluster ID
  status       Test status (tested / skipped)
  loglr        Log likelihood ratio
  df           Degrees of freedom
  p            Nominal p-value
  p.adjust     FDR-adjusted p-value (Benjamini-Hochberg)

--------------------------------------------------------------------------------
CD4Naive_sex_effect_sizes.txt
--------------------------------------------------------------------------------
Format   : Tab-delimited text, one row per intron
Content  : Effect sizes (delta PSI) for each intron in tested clusters

Columns  :
  intron       Intron ID
  logef        Log effect size
  male         Mean PSI in male group
  female       Mean PSI in female group
