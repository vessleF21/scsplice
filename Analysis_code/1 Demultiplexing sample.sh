


# Genotyping data processing: this script covers the pipeline for genotype generation, quality control, imputation, and post-imputation filtering

# Official links:
#   array-analysis-cli    : https://support.illumina.com/array/array_software/ima-array-analysis-cli.html
#   PLINK v1.9            : https://www.cog-genomics.org/plink/
#   bcftools              : https://samtools.github.io/bcftools/bcftools.html
#   bgzip / tabix (htslib): https://github.com/samtools/htslib


# =============================================================================
# Part A. Genotype Calling: converts raw IDAT intensity files to called genotypes
# =============================================================================

# A1: GenCall Text/Call（GTC）from IDAT files
# Reads raw fluorescence intensity from IDAT files and calls genotypes using the
# cluster file (EGT). Output is one GTC file per sample containing genotype calls
# and confidence scores. BPM manifest defines the probe sequences on the array.
# IDAT_DIR: directory containing raw IDAT files from the Illumina SNP array
# BPM: Bead Pool Manifest defining probe sequences and chromosome positions
# EGT: cluster file containing intensity cluster boundaries for genotype calling
$CLI genotype call \
     --bpm-manifest $BPM \
     --cluster-file $EGT \
     --idat-folder $IDAT_DIR \
     --output-folder gtc

# A2: GTC to called genotypes
# Converts each GTC file to a called genotypes by resolving alleles against the GRCh38
# reference genome. The CSV manifest provides strand and sequence context required
# for correct REF/ALT assignment. 
# CSV: CSV-format manifest containing SourceSeq required for allele resolution
# FASTA: GRCh38 reference genome; a .fai index must exist in the same directory
$CLI genotype gtc-to-vcf \
     --bpm-manifest $BPM \
     --csv-manifest $CSV \
     --genome-fasta-file $FASTA \
     --gtc-folder gtc \
     --output-folder geno

# =============================================================================
# Part B. SNP-level QC 
# SNPs are filtered in four sequential steps. Each step outputs a new PLINK
# binary dataset; wc -l on the .bim file reports how many SNPs remain after
# each filter
# =============================================================================

# B1. Restrict to autosomes (chromosomes 1-22)
plink \
    --bfile geno \
    --autosome \
    --make-bed \
    --out auto
wc -l geno.bim auto.bim

# B2. Remove SNPs with call rate < 95%
plink \
    --bfile auto \
    --geno 0.05 \
    --make-bed \
    --out snp_geno
wc -l auto.bim snp_geno.bim

# B3. Remove SNPs with minor allele frequency < 1%
plink \
    --bfile snp_geno \
    --maf 0.01 \
    --make-bed \
    --out snp_maf
wc -l snp_geno.bim snp_maf.bim

# B4. Remove SNPs deviating from Hardy-Weinberg equilibrium (P < 1e-6)
plink \
    --bfile snp_maf \
    --hwe 1e-6 \
    --make-bed \
    --out snpqc
wc -l snp_maf.bim snpqc.bim


# =============================================================================
# Part C. Sample-level QC
# Samples are filtered in four sequential steps targeting low-quality genotyping,
# abnormal heterozygosity, cryptic relatedness, and population stratification.
# =============================================================================

# C1. Remove samples with call rate < 95%
plink \
    --bfile snpqc \
    --mind 0.05 \
    --make-bed \
    --out mind

# C2. Remove heterozygosity outliers
plink \
    --bfile mind \
    --het \
    --out het

plink \
    --bfile mind \
    --remove het.outliers \
    --make-bed \
    --out het

# C3. Remove related individuals; retain higher call-rate sample per pair
plink \
    --bfile het \
    --indep-pairwise 50 5 0.2 \
    --exclude LD-regions.txt \
    --range \
    --out prune

plink \
    --bfile het \
    --extract prune.prune.in \
    --genome \
    --out ibd

plink \
    --bfile het \
    --remove related.remove \
    --make-bed \
    --out unrel

# C4. Remove ancestry outliers via PCA
plink \
    --bfile unrel \
    --extract prune.prune.in \
    --pca \
    --out pca

plink \
    --bfile unrel \
    --remove pca.outliers \
    --make-bed \
    --out clean

# =============================================================================
# Part D. Imputation
# Split clean genotypes by chromosome into bgzip-compressed VCF
# Each chromosome is exported as an individual VCF file and compressed with bgzip.
# =============================================================================

for chr in $(seq 1 22); do
  plink \
      --bfile clean \
      --chr $chr \
      --recode vcf-iid \
      --out chr$chr
  bcftools sort chr$chr.vcf | bgzip -c > chr$chr.vcf.gz
done

# =============================================================================
# Part E. Post-imputation filtering
# Filters: imputation quality R² > 0.8 | biallelic SNPs only | MAF >= 5% | HWE P >= 1e-6
# =============================================================================

# E1. Per-chromosome: keep R²>0.8, biallelic SNPs only (-m2 -M2 -v snps)
for chr in $(seq 1 22); do
  bcftools view \
      -i 'INFO/R2>0.8' \
      -m2 -M2 -v snps \
      chr$chr.dose.vcf.gz \
      -Oz \
      -o chr$chr.f.vcf.gz
  tabix -p vcf chr$chr.f.vcf.gz
done

bcftools concat \
    $(for c in $(seq 1 22); do echo chr$c.f.vcf.gz; done) \
    -Oz \
    -o merged.vcf.gz

# E2. Convert to PLINK, autosomes only
plink \
    --vcf merged.vcf.gz \
    --double-id \
    --autosome \
    --make-bed \
    --out post_imp
wc -l post_imp.bim

# E3. Remove SNPs with MAF < 5%
plink \
    --bfile post_imp \
    --maf 0.05 \
    --make-bed \
    --out post_maf
wc -l post_imp.bim post_maf.bim

# E4. Remove SNPs deviating from HWE (P < 1e-6)
plink \
    --bfile post_maf \
    --hwe 1e-6 \
    --make-bed \
    --out final

wc -l post_maf.bim final.bim



# =============================================================================
# scRNA-seq preprocessing: SRA to FASTQ conversion and Cell Ranger alignment
# This script covers two steps:
#   1. fasterq-dump : converts downloaded .sra files to FASTQ format
#   2. Cell Ranger  : aligns reads and generates expression matrix and BAM
#
# Official links:
#   SRA Toolkit  : https://github.com/ncbi/sra-tools
#   Cell Ranger  : https://www.10xgenomics.com/support/software/cell-ranger
# ----------------------------------------------------------------------------

# =============================================================================
# Part F. SRA -> FASTQ (fasterq-dump)
# Converts each .sra file to FASTQ format. 
#--include-technical retains barcode and UMI reads required for 10x Chromium data processing. 
#--split-files
# produces three files per run:
#   _1.fastq
#   _2.fastq
#   _3.fastq
# =============================================================================

for SRR in $(ls ${SRA_DIR}/*.sra | xargs -n1 basename | sed 's/\.sra//'); do
    echo "Processing ${SRR}..."
    $FASTERQ ${SRA_DIR}/${SRR}.sra \
        -O ${FASTQ_DIR} \
        --include-technical \
        --split-files \
        -e ${THREADS}
done

# =============================================================================
# Part G. FASTQ -> aligned BAM + expression matrix (Cell Ranger count)
# Cell Ranger aligns reads to GRCh38, corrects cell barcodes, deduplicates
# UMIs, and generates a feature-barcode matrix and coordinate-sorted BAM.
# The BAM file (possorted_genome_bam.bam) is required for downstream
# genotype-based demultiplexing (Demuxlet).
# =============================================================================
$CELLRANGER count \
    --id=$SAMPLE_ID \
    --transcriptome=$TRANSCRIPTOME \
    --fastqs=$FASTQ_DIR \
    --sample=$SAMPLES \
    --create-bam true \
    --nosecondary
# --id             : output folder name; results written to ./<SAMPLE_ID>/outs/
# --transcriptome  : path to Cell Ranger-compatible GRCh38 reference package
# --fastqs         : directory containing FASTQ files for all runs in this batch
# --sample         : comma-separated SRR run IDs; Cell Ranger matches by filename prefix
# --create-bam     : generate coordinate-sorted BAM (possorted_genome_bam.bam),
#                    required for Demuxlet downstream analyses
# --nosecondary    : skip Cell Ranger's own dimensionality reduction and clustering,
#                    as custom QC and analysis will be performed downstream


## Use a single batch of data to illustrate the complete demultiplexing workflow ##

# =============================================================================
# Part H. Sample Demultiplexing with Demuxlet
# Assigns droplets to individual samples using genotype posterior probability from exonic SNPs
# Classifies droplets as singlets, doublets, or ambiguous based on SNP coverage
# possorted_genome_bam.bam is the Cell Ranger-generated BAM file for a single batch of single-cell fastq data, for example GSM5899873, and impute_samples.vcf is the genotype VCF corresponding to that same batch.
# Official Popscle repository link: https://github.com/statgen/popscle
# =============================================================================

popscle demuxlet --sam possorted_genome_bam.bam \
                 --vcf impute_samples.vcf \
                 --field GT \
                 --out impute_te1

