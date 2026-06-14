================================================================================
cis-sQTL Pipeline - Input/Output Description
================================================================================

================================================================================
INPUT FILES
================================================================================

--------------------------------------------------------------------------------
1. chr1.dose.vcf.gz
--------------------------------------------------------------------------------
Format   : VCF (bgzip compressed)
Content  : Imputed genotype dosage values (0~2) per sample per SNP
Scale    : chr1, 98 samples

--------------------------------------------------------------------------------
2. PC.txt
--------------------------------------------------------------------------------
Format   : Tab-delimited text
Content  : Covariate matrix (genotype PCs, cell-type proportions, etc.)
Scale    : 98 samples
Columns  :
  Col 1   Covariate name (e.g. PC1, PC2, ...)
  Col 2+  Covariate value per sample

--------------------------------------------------------------------------------
3. phenotype.txt.gz
--------------------------------------------------------------------------------
Format   : BED (bgzip compressed)
Content  : Splicing junction ratios (PSI) per sample per junction
Scale    : chr1, 98 samples
Columns  :
  Col 1   #Chr      Chromosome
  Col 2   start     Start coordinate (0-based)
  Col 3   end       End coordinate
  Col 4   pid       Junction ID
  Col 5   gid       Gene/cluster ID
  Col 6   strand    Strand
  Col 7+  [sample]  Splicing ratio value per sample

================================================================================
OUTPUT FILES
================================================================================

--------------------------------------------------------------------------------
permutation.txt
--------------------------------------------------------------------------------
Content  : Raw permutation results, one line per gene
           Includes beta distribution parameters and empirical p-values

--------------------------------------------------------------------------------
permutations_full.txt.thresholds.txt
--------------------------------------------------------------------------------
Content  : Per-gene nominal p-value significance thresholds

--------------------------------------------------------------------------------
permutations_full.txt.significant.txt
--------------------------------------------------------------------------------
Content  : Genes passing FDR < 0.05

--------------------------------------------------------------------------------
conditions.txt
--------------------------------------------------------------------------------
Content  : Significant cis-sQTLs after conditional analysis
           One row per independent signal per gene

--------------------------------------------------------------------------------
nominals_100.txt
--------------------------------------------------------------------------------
Content  : First 100 lines of nominal pass output (head -n 100)
           Full nominal output truncated for file size
