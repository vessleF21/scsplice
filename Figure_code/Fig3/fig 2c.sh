
samtools merge onek1k.bam BMem.bam BNavie.bam CD4EM.bam CD4Naive.bam CD4Treg.bam CD8GZMH.bam CD8GZMK.bam CD8Naive.bam MAIT.bam cM.bam NKBright.bam NKDim.bam ncM.bam

samtools sort onek1k.bam -o onek1k.sort.bam
samtools index onek1k.sort.bam

bamCoverage -b onek1k.sort.bam -o onek1k.sort.bw --normalizeUsing RPKM -bs 5 -p 16

computeMatrix scale-regions \
              -R gencode.v32.gtf \
              -S OneK1K.sort.bw \
              --metagene \
              -o oneK1K.mat.gz --outFileSortedRegions sortedRegions.bed \
              -b 1000 -a 1000 \
              --sortUsing mean \
              --skipZeros \
              -p 16

plotHeatmap -m onek1k.mat.gz \
    --samplesLabel Read2 \
    --colorMap viridis \
    --sortUsing mean \
    --sortRegions descend \
    --plotTitle "OneK1K Metagene Heatmap" \
    --heatmapWidth 6 --heatmapHeight 9 \
    -out onek1k.pdf
