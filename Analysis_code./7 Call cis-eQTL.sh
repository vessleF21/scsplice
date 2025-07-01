


import sys
import rpy2.robjects as ro
import tensorqtl
import pandas as pd
import torch
import rpy2
from tensorqtl import pgen, cis, trans, post
from pyplink import PyPlink

celltype = sys.argv[1]

phenotype_bed_file = f"{celltype}_sorted.bed.gz"
covariates_file = f"{celltype}_ok_PC.txt"
pgen_file = f"{celltype}_tsst"

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

pgr = pgen.PgenReader(pgen_file)
genotype_df = pgr.load_genotypes()
variant_df = pgr.variant_df

variant_df['chrom'] = variant_df['chrom'].apply(lambda x: 'chr' + str(x))

cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df)
tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0.85)
post.calculate_qvalues(cis_df, fdr=0.05, qvalue_lambda=0.85)

cis_df.to_csv(f'{celltype}_cis_df_saved.txt', sep='\t', index=True)

prefix = f"{celltype}_cisqtl_sample"
cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix, covariates_df, output_dir='.')

parquet_files = [f"{prefix}.cis_qtl_pairs.chr{i}.parquet" for i in range(1, 23)]
for parquet_file in parquet_files:
    try:
        data = pd.read_parquet(parquet_file)
        txt_file = parquet_file.replace('.parquet', '.txt')
        data.to_csv(txt_file, sep='\t', index=False)
        print(f"Converted {parquet_file} to {txt_file}")
    except Exception as e:
        print(f"Failed to process {parquet_file}: {e}")

indep_df = cis.map_independent(genotype_df, variant_df, cis_df, phenotype_df, phenotype_pos_df, covariates_df)
indep_df.to_csv(f'{celltype}_indep_df_saved.txt', sep='\t', index=True)

trans_df = trans.map_trans(genotype_df, phenotype_df, covariates_df, return_sparse=True, pval_threshold=1e-5, maf_threshold=0.05, batch_size=20000)
trans_df = trans.filter_cis(trans_df, phenotype_pos_df, variant_df, window=5000000)
trans_df.to_csv(f"{celltype}_trans_df.csv", index=False)










