================================================================================
sQTL-GWAS Colocalization Analysis - Input/Output Description
================================================================================

================================================================================
OVERVIEW
================================================================================
Tests for shared causal variants between sQTLs and GWAS traits using the
COLOC approximate Bayes factor framework. For each intron cluster, computes
posterior probabilities for five colocalization hypotheses (H0-H4).

================================================================================
DIRECTORY STRUCTURE
================================================================================

├── gwas_1/
│   └── CD4Naive_chr1_nominals_1/
│       ├── LiuZ_NatGenet_2023_CD/                <- case-control GWAS (Crohn's disease)
│       │   └── chr1_1000000_1005000_clu_1
│       └── ChenMH_CELL_2020_LYM/                 <- quantitative GWAS (Lymphocyte count)
│           └── chr1_1000000_1005000_clu_1
├── sqtl_1/
│   └── CD4Naive_chr1_nominals_1/
│       ├── LiuZ_NatGenet_2023_CD/
│       │   └── chr1_1000000_1005000_clu_1
│       └── ChenMH_CELL_2020_LYM/
│           └── chr1_1000000_1005000_clu_1
└── MAF/
    └── split_files/
        └── chrchr1.txt

================================================================================
INPUT FILES
================================================================================

--------------------------------------------------------------------------------
1. GWAS summary statistics
   gwas_1/{celltype}_{chr}_nominals_1/{gwas_name}/{intron}
   e.g. gwas_1/CD4Naive_chr1_nominals_1/LiuZ_NatGenet_2023_CD/chr1_1000000_1005000_clu_1
--------------------------------------------------------------------------------
Format   : Tab-separated, with header
Content  : GWAS summary statistics for one intron locus

Columns  :
  Col 1   dummy index (ignored)
  Col 2   variant_id      SNP ID (chr:pos:ref:alt)
  Col 3   beta            Effect size
  Col 4   se              Standard error
  Col 5   pval_nominal    Nominal p-value

Trait types supported:
  Case-control (cc)  : e.g. LiuZ_NatGenet_2023_CD (Crohn's disease)
  Quantitative       : e.g. ChenMH_CELL_2020_LYM  (Lymphocyte count)

--------------------------------------------------------------------------------
2. sQTL summary statistics
   sqtl_1/{celltype}_{chr}_nominals_1/{gwas_name}/{intron}
   e.g. sqtl_1/CD4Naive_chr1_nominals_1/LiuZ_NatGenet_2023_CD/chr1_1000000_1005000_clu_1
--------------------------------------------------------------------------------
Format   : Space-separated, no header
Content  : sQTL nominal pass results for one intron cluster

Columns  :
  Col 1   gene            Ensembl gene ID
  Col 2   intron          Intron cluster ID
  Col 3   variant_id      SNP ID (chr:pos:ref:alt)
  Col 4   pval_nominal    Nominal p-value
  Col 5   beta            Effect size (PSI change per allele)

--------------------------------------------------------------------------------
3. MAF file
   MAF/split_files/chr{chr}.txt
   e.g. MAF/split_files/chrchr1.txt
--------------------------------------------------------------------------------
Format   : Tab-separated, with header
Content  : Minor allele frequencies for all SNPs on the chromosome

Columns  :
  Col 1   variant_id      SNP ID (must match GWAS and sQTL files)
  Col 2   maf             Minor allele frequency (0.05 - 0.50)

================================================================================
OUTPUT FILES
================================================================================

--------------------------------------------------------------------------------
{celltype}_{chr}_nominals_1/{gwas_name}/{intron}
e.g. CD4Naive_chr1_nominals_1/LiuZ_NatGenet_2023_CD/chr1_1000000_1005000_clu_1
--------------------------------------------------------------------------------
Format   : Tab-separated, no header, 16 columns, one row per intron cluster
Content  : Colocalization results for this locus

Columns  :
  1   sqtl_variant_id_most    Top sQTL variant ID (before merging)
  2   intron_id               Intron cluster ID
  3   gene                    Gene ID
  4   nsnps                   Number of SNPs tested
  5   PP.H0.abf               Posterior prob: no association with either trait
  6   PP.H1.abf               Posterior prob: GWAS association only
  7   PP.H2.abf               Posterior prob: sQTL association only
  8   PP.H3.abf               Posterior prob: both traits, different causal variants
  9   PP.H4.abf               Posterior prob: both traits, shared causal variant
  10  sqtl_variant_P_most     Top sQTL p-value (before merging)
  11  sqtl_variant_id_merge   Top sQTL variant ID (after merging with GWAS)
  12  sqtl_variant_P_merge    Top sQTL p-value (after merging)
  13  GWAS_variant_id_most    Top GWAS variant ID (before merging)
  14  GWAS_variant_P_most     Top GWAS p-value (before merging)
  15  GWAS_variant_id_merge   Top GWAS variant ID (after merging)
  16  GWAS_variant_P_merge    Top GWAS p-value (after merging)

--------------------------------------------------------------------------------
coloc_results.txt
--------------------------------------------------------------------------------
Format   : Plain text
Content  : Terminal output from both COLOC runs (PP values printed to stdout)
