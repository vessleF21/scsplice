


computeMatrix scale-regions \
              -R gencode.v32.gtf \
              -S OneK1K.sort.bw \
              --metagene \
              -o oneK1K.mat.gz --outFileSortedRegions sortedRegions.bed \
              -b 1000 -a 1000 \
              --sortUsing mean \
              --skipZeros \
              -p 16


plotHeatmap -m onek1k001.final.mat.gz \
    --samplesLabel Read2 \
    --colorMap viridis \
    --sortUsing mean \
    --sortRegions descend \
    --plotTitle "OneK1K Metagene Heatmap" \
    --heatmapWidth 6 --heatmapHeight 9 \
    -out onek1k001.pdf