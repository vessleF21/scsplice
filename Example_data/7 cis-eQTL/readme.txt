================================================================================
cis-eQTL Mapping Pipeline - Input/Output Description
================================================================================

================================================================================
INPUT FILES
================================================================================

--------------------------------------------------------------------------------
MAIT_sorted_chr22_200.bed.gz
--------------------------------------------------------------------------------
Format   : BED (bgzip compressed, tabix indexed)
Content  : Pseudobulk gene expression matrix after normalization
           Subset to chr22 and 200 samples for testing

Columns  :
  Col 1   #chr        Chromosome (chr22)
  Col 2   start       Gene start position
  Col 3   end         Gene end position
  Col 4   gene_id     Ensembl gene ID (version number removed)
  Col 5+  [sample]    Inverse-normal transformed expression per donor

Scale    : chr22 only; 200 donors; genes detected in >= 40% of samples

--------------------------------------------------------------------------------
MAIT_ok_PC_200.txt
--------------------------------------------------------------------------------
Format   : Tab-delimited text
Content  : Covariate matrix combining PEER factors and sample metadata
           Subset to 200 samples matching expression BED

Rows     : Covariates (PEER factors + genotype PCs + sex + age)
Columns  : Donors (must match column order of expression BED)

  InferredCov1..N    PEER factors (N depends on sample size)
  PC1..7             Genotype principal components
  sex                Biological sex (encoded as numeric)
  age                Age at sample collection

--------------------------------------------------------------------------------
MAIT_chr22_200.pgen
MAIT_chr22_200.pvar
MAIT_chr22_200.psam
--------------------------------------------------------------------------------
Format   : PLINK2 binary genotype format

  MAIT_chr22_200.pgen    Genotype matrix (binary)
  MAIT_chr22_200.pvar    Variant information (chr, pos, ref, alt, ID)
  MAIT_chr22_200.psam    Sample information

Scale    : chr22 only; 200 donors matching expression BED and PC file

================================================================================
OUTPUT FILES
================================================================================

--------------------------------------------------------------------------------
MAIT_chr22_cis_df_saved.txt
--------------------------------------------------------------------------------
Format   : Tab-delimited text
Content  : cis permutation pass results, one row per gene
           Includes empirical p-values and FDR q-values

Columns  :
  phenotype_id       Gene ID
  variant_id         Lead variant ID
  tss_distance       Distance from TSS to lead variant
  beta               Effect size (normalized)
  pval_nominal       Nominal p-value of lead variant
  pval_perm          Permutation p-value
  pval_beta          Beta-approximated p-value
  qval               FDR q-value (qvalue_lambda=0.85)

--------------------------------------------------------------------------------
MAIT_chr22_cisqtl_sample.cis_qtl_pairs.chr22.parquet
--------------------------------------------------------------------------------
Format   : Apache Parquet
Content  : cis nominal pass results; all tested variant-gene pairs on chr22

Columns  :
  phenotype_id       Gene ID
  variant_id         Variant ID
  tss_distance       Distance from TSS
  maf                Minor allele frequency
  beta               Effect size
  se                 Standard error
  pval_nominal       Nominal p-value

--------------------------------------------------------------------------------
MAIT_chr22_indep_df_saved.txt
--------------------------------------------------------------------------------
Format   : Tab-delimited text
Content  : Conditionally independent cis-eQTLs per gene
           Identifies secondary signals after conditioning on lead variant

Columns  :
  phenotype_id       Gene ID
  variant_id         Independent variant ID
  rank               Signal rank (0=lead, 1=secondary, ...)
  beta               Effect size
  pval_nominal       Nominal p-value
  pval_beta          Beta-approximated p-value
