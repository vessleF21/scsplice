

rm(list=ls())
gc()

library(Seurat)
library(patchwork)
library(harmony)
library(tidyverse)

options(future.globals.maxSize = 1000 * 1024^3)

source("custom_seurat_functions.R")

sce = readRDS("T_sce.qs")

colorrSS = c("#4E79A7","#E15759","#499894","#F28E2B","#8CD17D","#FFBE7D","#B6992D","#FF9D9A","#D37295","#A0CBE8","#79706E", "#5D69B1", "#24796C","#499894",'#DAA51B', '#000000', '#99C945', '#ED645A',"#7B6FD0", "#CF4A31", "#D0CD47","#722A2D", "#CBC594", "#D19EC4")
colorrSS = colorrSS[1:length(levels(Idents(sce)))]
names(colorrSS) = levels(Idents(sce))

p1 <- DimPlot(sce, reduction = "umap",
              group.by = "test", 
              label = T,label.box = T,repel = T,raster = F)&ggtitle("T_NK")&
  scale_color_manual(values = colorrSS)&
  scale_fill_manual(values = colorrSS)

ggsave(plot= p1, filename = "Annotation_Umap_Vis.png", 
       width = 5.2, height = 4.5)

check_genes = c('CD3D','CD3E',  #T cells
                'CD4','CD40LG', #CD4T cells
                "LEF1","CCR7", "TCF7", # Naive
                "TNFRSF4","GPR183","TIMP1", #Tem
                'FOXP3',"RTKN2",'IL2RA', #Treg
                'CD8A','CD8B', #CD8T
                "GZMH", # GZMH+ CD8T
                "GZMK","CXCR4","CCL5", # GZMK+ CD8T
                "SLC4A10","ZBTB16","RORC", #C3_CD8-SLC4A10 mucosal-associated invariant T cells (MAIT)，MAIT
                "FCGR3A","NCAM1","KLRC1",   # NK1 CD16_neg bright, FCGR3A low NCAM1 bright
                "TRDC" #progenitor cells (Progen) in PBMC PMID: 35389781
)

sce$test = factor(sce$test, levels = c("CD4Naive", "CD4EM", "CD4Treg",'CD8Naive', "CD8GZMH", "CD8GZMK", "MAIT", "NKDim", "NKBright","gdT"))
Idents(sce) = "test"


p2 = DotPlot_2(object = sce,
               features = check_genes,
               assay = "RNA", scale = T,
               Combine = F, color.use = colorrSS)&
  scale_color_distiller(palette = 'RdYlBu')

ggsave(plot= p2, filename = "Annotation_Markers_Vis.pdf", 
       width = 12, height = 6, units = "cm")










