


rm(list=ls())
gc()

library(Seurat)
library(patchwork)
library(harmony)
library(tidyverse)

options(future.globals.maxSize = 1000 * 1024^3)

source("custom_seurat_functions.R")

sce = readRDS("B_sce.qs")
colorrSS = c("#4E79A7","#E15759","#499894","#F28E2B","#8CD17D","#FFBE7D","#B6992D","#FF9D9A","#D37295","#A0CBE8","#79706E")
colorrSS = colorrSS[1:length(levels(Idents(sce)))]
names(colorrSS) = levels(Idents(sce))

p1 <- DimPlot(sce, reduction = "umap",
              group.by = "CT_res", 
              label = T,label.box = T,repel = T,raster = F)&ggtitle("B")&
  scale_color_manual(values = colorrSS)&
  scale_fill_manual(values = colorrSS)

ggsave(plot= p1, filename = "Annotation_Umap_Vis.png", 
       width = 5.2, height = 4.5)


check_genes = c("IL4R","FCER2","TCL1A", #naive B cells
                "TNFRSF13B","AIM2", "BANK1", #Bmem
                "FCRL5", "TBX21", "CD19", "MS4A1", #BAtypical
                "MZB1","JCHAIN" #B plasam 
)
p2 = DotPlot_2(object = sce,
               features = check_genes,
               assay = "RNA", scale = T,
               Combine = F, color.use = colorrSS)&
  scale_color_distiller(palette = 'RdYlBu')

ggsave(plot= p2, filename = "Annotation_Markers_Vis.pdf", 
       width = 12, height = 6, units = "cm")