
/public/home/liangyan/software/QTLtools trans \
      --vcf chr1_22.0.9.vcf.gz \
      --bed phenotype.txt.gz \
      --cov PC.txt \
      --out trans.sample \
      --sample 5000 \
      --threshold 1e-5 --normal --window 5000000 --bin 1000

/public/home/liangyan/software/QTLtools trans \
          --vcf chr1_22.0.9.vcf.gz \
          --bed phenotype.txt.gz \
          --cov PC.txt \
          --out trans.adjust \
          --adjust trans.best.txt.gz \
          --threshold 0.1 --normal --window 5000000 --bin 1000

Rscript runFDR.mappability.R celltype_chr_trans.txt.gz \
            0.05 trans.adjust.txt \
           output


           