
library(optparse)
option_list = list(
  make_option("--filename", action="store", default=NA, type='character', help="Path to filename [required]"),
  make_option("--batchss", action="store", default=NA, type='character', help="batch for umap [required]"),
  make_option("--dimss", action="store", default=0, type='integer', help="dimss for number"),
  make_option("--cellTY", action="store", default=NA, type='character', help="batch for umap [required]"),
  make_option("--filter", action="store_true", default=FALSE, help="filter or not"),
  make_option("--neighbors_within_batch", action="store_true", default=FALSE, help="neighbors_within_batch or not")
)

opt = parse_args(OptionParser(option_list=option_list))

library(Seurat)
library(patchwork)
library(harmony)
library(tidyverse)
library(clustree)
library(cowplot)
library(qs)
library(bbknnR)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
options(future.globals.maxSize = 1000 * 1024^3)
sce = qread(file = opt$filename)
index <- match(paste0("RNA_snn_res.", seq(0.1, 3, by=0.1)), colnames(sce@meta.data))
index <- index[!is.na(index)]
sce@meta.data <- sce@meta.data[,-index]

#sce@meta.data$ct_cov[sce@meta.data$ct_cov == ""] <- "unknown"

sce <- NormalizeData(object = sce, normalization.method = "LogNormalize", scale.factor = 1e4)

sce <- FindVariableFeatures(sce)

sce <- ScaleData(object = sce, features = rownames(sce))
sce <- RunPCA(object = sce, features = VariableFeatures(object = sce), npcs = 50)

sce <- sce %>% 
  RunHarmony(group.by.vars = "POOL")

sce <- sce %>% 
  FindNeighbors(reduction = "harmony", dims = 1:opt$dimss)%>% 
  RunUMAP(reduction = "harmony", dims = 1:opt$dimss, verbose = F)


p1 <- DimPlot(sce, reduction = "umap", group.by = "POOL", raster = F)
p2 <- DimPlot(sce, reduction = "umap",group.by = "l2_onek1k", label = TRUE, raster = F)
p3 <- DimPlot(sce, reduction = "umap",group.by = "cell_type_onek1k", label = TRUE, raster = F)
P.total = wrap_plots(p1,p2,p3)
ggsave(plot = P.total, filename = paste0(opt$cellTY,"_",opt$dimss,"_",opt$batchss,"_Step6.Umap_plot.png"),width = 48,height = 6)

qsave(sce, file = paste0(opt$cellTY,"_",opt$dimss,"_",opt$batchss,"_before_sce.harmony.qs"))

for (res in c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)) {
  print(res)
  sce <- FindClusters(sce,resolution = res, 
                      algorithm = 1)
}
qsave(sce, file=paste0(opt$cellTY,"_",opt$dimss,"_",opt$batchss,"_batch_cov_sce_qc_after_findcluster.qs"))

cluster_umap <- wrap_plots(ncol = 4,
                           DimPlot(sce, reduction = "umap", 
                                   group.by = "RNA_snn_res.0.3", label = T,raster=FALSE) & NoAxes(),
                           DimPlot(sce, reduction = "umap", 
                                   group.by = "RNA_snn_res.0.4", label = T,raster=FALSE) & NoAxes(),
                           DimPlot(sce, reduction = "umap", 
                                   group.by = "RNA_snn_res.0.5", label = T,raster=FALSE) & NoAxes(), 
                           DimPlot(sce, reduction = "umap", 
                                   group.by = "RNA_snn_res.0.6", label = T,raster=FALSE) & NoAxes(),
                           DimPlot(sce, reduction = "umap", 
                                   group.by = "RNA_snn_res.0.7", label = T,raster=FALSE)& NoAxes(),
                           DimPlot(sce, reduction = "umap", 
                                   group.by = "RNA_snn_res.0.8", label = T,raster=FALSE) & NoAxes(),
                           DimPlot(sce, reduction = "umap", 
                                   group.by = "RNA_snn_res.0.9", label = T,raster=FALSE) & NoAxes(), 
                           DimPlot(sce, reduction = "umap", 
                                   group.by = "RNA_snn_res.1", label = T,raster=FALSE) & NoAxes(),
                           DimPlot(sce, reduction = "umap", 
                                   group.by = "RNA_snn_res.1.1", label = T,raster=FALSE) & NoAxes(), 
                           DimPlot(sce, reduction = "umap", 
                                   group.by = "RNA_snn_res.1.2", label = T,raster=FALSE) & NoAxes()
)
ggsave(cluster_umap,filename = paste0(opt$cellTY,"_",opt$dimss,"_",opt$batchss,"_batch_cov","_BBKNN_Umap_plot_2.png"),width = 30, height = 10)