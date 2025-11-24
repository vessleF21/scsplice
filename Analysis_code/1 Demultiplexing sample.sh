
## Use a single batch of data to illustrate the complete demultiplexing workflow ##

step 1: Sample Demultiplexing with Demuxlet
# Assigns droplets to individual samples using genotype posterior probability from exonic SNPs
# Classifies droplets as singlets, doublets, or ambiguous based on SNP coverage
# possorted_genome_bam.bam is the Cell Ranger–generated BAM file for a single batch of single-cell fastq data, for example GSM5899873, and impute_samples.vcf is the genotype VCF corresponding to that same batch.
# Official Popscle repository link: https://github.com/statgen/popscle

popscle demuxlet --sam possorted_genome_bam.bam \
                 --vcf impute_samples.vcf \
                 --field GT \
                 --out impute_te1



step 2: BAM File Splitting by Sample
# Filter and split BAM file by sample-specific barcodes
# Extracts reads belonging to singlet cells assigned to single sample-onek1k_sample1
# Uses 20 threads for parallel processing
# Official sinto repository link: https://github.com/timoast/sinto
sinto filterbarcodes -p 20 \
                     -b possorted_genome_bam.bam \
                     -c onek1k_sample1.txt \
                     --outdir ./split_bam/onek1k_sample1


step 3: VCF File Processing for Individual Samples
# Extract sample-specific variants with at least 1 alternate allele
# Removes monomorphic sites for the specific sample
# Official bcftools repository link: https://github.com/samtools/bcftools
bcftools view -c1 -s onek1k_sample1 onek1k_sample1_vcf

# Sort VCF file and compress with bgzip
bcftools sort -Oz -o onek1k_sample1.vcf.gz

# Create index for compressed VCF
bcftools index onek1k_sample1.vcf.gz


step 4: STAR Alignment with WASP Filtering
# Two-pass STAR alignment with single-cell and allele-specific mapping features
# Key parameters:
# - Uses GRCh38 reference genome
# - Forward strand specificity for 10x Chromium chemistry
# - CB_UMI_Simple mode for 10x data processing
# - 1MM_Directional_UMItools for UMI deduplication
# - WASP filtering to remove reference allele mapping bias
# - Outputs SAM attributes required for downstream QC and quantification
# Official STAR repository link: https://github.com/alexdobin/STAR

STAR \
    --runThreadN 24 \
    --genomeDir ./star/hg38_gencode_v32 \
    --twopassMode Basic \
    --soloStrand Forward \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist ./ref/737K-august-2016.txt \
    --soloCBmatchWLtype 1MM \
    --soloUMIdedup 1MM_Directional_UMItools \
    --readFilesIn $readFilesIn \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM XS vG vA vW \
    --waspOutputMode SAMtag \
    --varVCFfile <(zcat onek1k_sample1.vcf.gz) \
    --sjdbGTFfile ./hg38/gencode.annotation.gtf \
    --soloOutFileNames onek1k_sample1 \
    --outFileNamePrefix onek1k_sample1


step 5: PCR Duplicate Removal
# BARCODE_TAG=CB: Uses cell barcode to identify duplicates within each cell
# MOLECULAR_IDENTIFIER_TAG=UB: Uses UMI-corrected barcode for molecular duplicate detection
# REMOVE_DUPLICATES=true: Physically removes duplicates from output BAM
# Official STAR repository link: https://broadinstitute.github.io/picard/

java -jar picard.jar MarkDuplicates \
    REMOVE_DUPLICATES=true \
    BARCODE_TAG=CB \
    MOLECULAR_IDENTIFIER_TAG=UB \
    I=onek1k_sample1.sorted.bam \
    O=onek1k_sample1.sorted_rmdup.bam \
    M=marked_dup_metrics.onek1k_sample1.txt

