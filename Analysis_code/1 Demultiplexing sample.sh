popscle  demuxlet --sam possorted_genome_bam.bam --vcf impute_samples.vcf --field GT --out impute_te1

popscle dsc-pileup --sam possorted_genome_bam.bam --vcf impute_samples.vcf --out impute_te2
popscle demuxlet --plp impute_te2 --vcf impute_samples.vcf --field GT --out impute_te2

sinto filterbarcodes -p 20 -b possorted_genome_bam.bam -c GSM5899873.txt --outdir ./split_bam/GSM5899873

bcftools view -c1 -s GSM5899873 GSM5899873_vcf 
bcftools sort -Oz -o GSM5899873.vcf.gz
bcftools index GSM5899873.vcf.gz

STAR \
        --runThreadN 24 \
        --genomeDir ./star/hg38_gencode_v32 \
        --twopassMode Basic \
        --soloStrand Forward \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist ./ref/737K-august-2016.txt --soloCBmatchWLtype 1MM --soloUMIdedup 1MM_Directional_UMItools \
        --readFilesIn $readFilesIn \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM XS vG vA vW \
        --waspOutputMode SAMtag --varVCFfile <(zcat GSM5899873.vcf.gz) \
        --sjdbGTFfile ./hg38/gencode.annotation.gtf \
        --soloOutFileNames GSM5899873 \
        --outFileNamePrefix GSM5899873

java -jar picard.jar MarkDuplicates \
          REMOVE_DUPLICATES=true \
          BARCODE_TAG=CB \
          MOLECULAR_IDENTIFIER_TAG=UB \
          I=GSM5899873.sorted.bam \
          O=GSM5899873.sorted_rmdup.bam \
          M=marked_dup_metrics.GSM5899873.txt