
================================================================================
trans-sQTL Test Dataset - Input/Output File Description
================================================================================

================================================================================
INPUT FILES
================================================================================

--------------------------------------------------------------------------------
1. --vcf  chr1.dose.vcf.gz (Please download it from Zenodo: https://zenodo.org/records/17979485)
--------------------------------------------------------------------------------
Format   : VCF (bgzip compressed + .tbi index)
Content  : Imputed genotype dosage values (0~2) per sample per SNP
Scale    : chr1, ~380,000 variants, 800 samples
Columns  : #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  [sample1] ...

--------------------------------------------------------------------------------
2. --bed  ct_phenotype.txt.gz
--------------------------------------------------------------------------------
Format   : BED (bgzip compressed + .tbi index, sorted by coordinate)
Content  : Splicing junction ratios (PSI) per sample per junction
Scale    : chr1, 100 genes (randomly selected)
Columns  :
  Col 1   #Chr      Chromosome
  Col 2   start     Start coordinate (0-based)
  Col 3   end       End coordinate
  Col 4   pid       Junction ID (e.g. chr1:1000:2000:clu_1_-)
  Col 5   gid       Gene/cluster ID
  Col 6   strand    Strand
  Col 7+  [sample]  Splicing ratio value per sample

--------------------------------------------------------------------------------
3. --cov  ct_PC.txt
--------------------------------------------------------------------------------
Format   : Tab-delimited text
Content  : Covariate matrix for confounding factor correction
           (genotype PCs, cell-type proportions, etc.)
Scale    : 800 samples
Columns  :
  Col 1   Covariate name (e.g. PC1, PC2, ...)
  Col 2+  Covariate value per sample
Note     : Sample IDs must match VCF and BED exactly

================================================================================
OUTPUT FILES
================================================================================

--------------------------------------------------------------------------------
permutation
--------------------------------------------------------------------------------
ct_trans.sample.best.txt.gz    Best trans hit per phenotype (null distribution)
ct_trans.sample.hits.txt.gz    All hits passing threshold
ct_trans.sample.bins.txt.gz    Bin statistics for FDR calibration

--------------------------------------------------------------------------------
condition
--------------------------------------------------------------------------------
ct_trans.adjust.best.txt.gz    Best trans hit per phenotype (adjusted p-values)
ct_trans.adjust.hits.txt.gz    All significant hits after adjustment
ct_trans.adjust.bins.txt.gz    Bin statistics

--------------------------------------------------------------------------------
FDR
--------------------------------------------------------------------------------
ct_trans.adjust.best.genes.txt          All genes with FDR-corrected p-values
ct_trans.adjust.best.genes.txt.filtered.txt   *** FINAL RESULT: significant trans-sQTLs ***

