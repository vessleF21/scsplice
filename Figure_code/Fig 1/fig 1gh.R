rm(list=ls())
gc()

library(Seurat)
library(patchwork)
library(harmony)
library(tidyverse)

options(future.globals.maxSize = 1000 * 1024^3)

source("custom_seurat_functions.R")

sce = readRDS("Mono_sce.qs")

colorrSS = c("#4E79A7","#E15759","#499894","#F28E2B","#8CD17D","#FFBE7D","#B6992D","#FF9D9A","#D37295","#A0CBE8","#79706E")
colorrSS = colorrSS[1:length(levels(Idents(sce)))]
names(colorrSS) = levels(Idents(sce))

p1 <- DimPlot(sce, reduction = "umap",
              group.by = "CT_res", 
              label = T,label.box = T,repel = T,raster = F)&ggtitle("Myeloid")&
  scale_color_manual(values = colorrSS)&
  scale_fill_manual(values = colorrSS)

ggsave(plot= p1, filename = "Annotation_Umap_Vis.png", 
       width = 5.2, height = 4.5)

check_genes = c("CD14", "S100A8","S100A9","FCN1", #CD14+ monocyte 
                "FCGR3A","CTSS","CDKN1C",#CD16+ monocyte 
                "C1QB","C1QC","C1QA",#ncM comp
                "CLEC9A","IDO1","IDO2", #cDC1 
                "FCER1A","CLEC10A","CD1C","CD1E", #cDC2
                "LILRA4","SHD","LRRC26","PACSIN1","IL3RA"#pDCs C9_pDC_LILRA4 GZMB+
)
Idents(sce) = "CT_res"
p2 = DotPlot_2(object = sce,
               features = check_genes,
               assay = "RNA", scale = T,
               Combine = F, color.use = colorrSS)&
  scale_color_distiller(palette = 'RdYlBu')

ggsave(plot= p2, filename = "Annotation_Markers_Vis.pdf", 
       width = 12, height = 6, units = "cm")
