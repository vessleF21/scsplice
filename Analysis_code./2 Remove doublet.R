
library(optparse)
option_list = list(
  make_option("--GSEname", action="store", default=NA, type='character', help="Path to filename [required]"),
  make_option("--Batchname", action="store", default=NA, type='character', help="Path to filename [required]"),
  make_option("--doublerate", action="store", default=NA, type='numeric', help="Path to filename [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(Seurat)
library(DoubletFinder)
library(harmony)
library(tidyverse)
library(data.table)
library(stringr)
options(future.globals.maxSize = 1000 * 1024^3)

source("doubletFinder.R")

paramSweep <- function(seu, PCs=1:10, sct = FALSE, num.cores=1) {
  require(Seurat); require(fields); require(parallel)
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05,0.3,by=0.05)

  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@meta.data)/(1-0.05) - nrow(seu@meta.data))
  pK.test <- round(pK*min.cells)
  pK <- pK[which(pK.test >= 1)]

  orig.commands <- seu@commands

  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000, replace=FALSE)]
    data <- seu@assays$RNA@counts[ , real.cells]
    n.real.cells <- ncol(data)
  }

  if (nrow(seu@meta.data) <= 10000){
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }

  if(num.cores>1){
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)),
                        FUN = parallel_paramSweep,
                        n.real.cells,
                        real.cells,
                        pK,
                        pN,
                        data,
                        orig.commands,
                        PCs,
                        sct,mc.cores=num.cores)
    stopCluster(cl)
  }else{
    output2 <- lapply(as.list(1:length(pN)),
                      FUN = parallel_paramSweep,
                      n.real.cells,
                      real.cells,
                      pK,
                      pN,
                      data,
                      orig.commands,
                      PCs,
                      sct)
  }

  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for(i in 1:length(output2)){
    for(j in 1:length(output2[[i]])){
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }

  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" ))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)

}


find.pK <- function(sweep.stats) {

  ## Implementation for data without ground-truth doublet classifications
  '%ni%' <- Negate('%in%')
  if ("AUC" %ni% colnames(sweep.stats) == TRUE) {
    ## Initialize data structure for results storage
    bc.mvn <- as.data.frame(matrix(0L, nrow=length(unique(sweep.stats$pK)), ncol=5))
    colnames(bc.mvn) <- c("ParamID","pK","MeanBC","VarBC","BCmetric")
    bc.mvn$pK <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)

    ## Compute bimodality coefficient mean, variance, and BCmvn across pN-pK sweep results
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind, "BCreal"])^2)
    }

    ## Plot for visual validation of BCmvn distribution
    #par(mar=rep(1,4))
    #x <- plot(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, pch=16, col="#41b6c4", cex=0.75)
    #x <- lines(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, col="#41b6c4")
    #print(x)

    return(bc.mvn)
  }
}

parallel_paramSweep <- function(n, n.real.cells, real.cells, pK, pN, data, orig.commands, PCs, sct)  {

  sweep.res.list = list()
  list.ind = 0

  ## Make merged real-artifical data
  print(paste("Creating artificial doublets for pN = ", pN[n]*100,"%",sep=""))
  n_doublets <- round(n.real.cells/(1 - pN[n]) - n.real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[, real.cells1] + data[, real.cells2])/2
  colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
  data_wdoublets <- cbind(data, doublets)

  ## Pre-process Seurat object
  if (sct == FALSE) {
    print("Creating Seurat object...")
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

    print("Normalizing Seurat object...")
    seu_wdoublets <- NormalizeData(seu_wdoublets,
                                   normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                   scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                   margin = orig.commands$NormalizeData.RNA@params$margin)

    print("Finding variable genes...")
    seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
                                          selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                          loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                          clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                          mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                          dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                          num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                          binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                          nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                          mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                          dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)

    print("Scaling data...")
    seu_wdoublets <- ScaleData(seu_wdoublets,
                               features = orig.commands$ScaleData.RNA$features,
                               model.use = orig.commands$ScaleData.RNA$model.use,
                               do.scale = orig.commands$ScaleData.RNA$do.scale,
                               do.center = orig.commands$ScaleData.RNA$do.center,
                               scale.max = orig.commands$ScaleData.RNA$scale.max,
                               block.size = orig.commands$ScaleData.RNA$block.size,
                               min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)

    print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets,
                            features = orig.commands$ScaleData.RNA$features,
                            npcs = length(PCs),
                            rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                            weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                            verbose=FALSE)
  }

  if (sct == TRUE) {
    require(sctransform)
    print("Creating Seurat object...")
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

    print("Running SCTransform...")
    seu_wdoublets <- SCTransform(seu_wdoublets)

    print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
  }

  ## Compute PC distance matrix
  print("Calculating PC distance matrix...")
  nCells <- nrow(seu_wdoublets@meta.data)
  pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
  rm(seu_wdoublets)
  gc()
  dist.mat <- fields::rdist(pca.coord)[,1:n.real.cells]

  ## Pre-order PC distance matrix prior to iterating across pK for pANN computations
  print("Defining neighborhoods...")
  for (i in 1:n.real.cells) {
    dist.mat[,i] <- order(dist.mat[,i])
  }

  ## Trim PC distance matrix for faster manipulations
  ind <- round(nCells * max(pK))+5
  dist.mat <- dist.mat[1:ind, ]

  ## Compute pANN across pK sweep
  print("Computing pANN across all pK...")
  for (k in 1:length(pK)) {
    print(paste("pK = ", pK[k], "...", sep = ""))
    pk.temp <- round(nCells * pK[k])
    pANN <- as.data.frame(matrix(0L, nrow = n.real.cells, ncol = 1))
    colnames(pANN) <- "pANN"
    rownames(pANN) <- real.cells
    list.ind <- list.ind + 1

    for (i in 1:n.real.cells) {
      neighbors <- dist.mat[2:(pk.temp + 1),i]
      pANN$pANN[i] <- length(which(neighbors > n.real.cells))/pk.temp
    }

    sweep.res.list[[list.ind]] <- pANN

  }

  return(sweep.res.list)
}


removeDoublet = function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs = 30)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  seurat_obj <- FindNeighbors(seurat_obj)
  seurat_obj <- FindClusters(seurat_obj)
  
  seurat_test <- paramSweep(seurat_obj, PCs = 1:30, sct = F)
  seurat_stats <- summarizeSweep(seurat_test, GT = FALSE)
  bcmvn <- find.pK(seurat_stats)
  
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  DoubletRate = as.numeric(opt$doublerate)
  homotypic.prop <- modelHomotypic(seurat_obj$seurat_clusters)
  nExp_poi <- round(DoubletRate*ncol(seurat_obj))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  
  seurat_obj <- doubletFinder(seurat_obj,
                                 pK = pK_bcmvn,
                                 PCs = 1:30,
                                 nExp = nExp_poi,
                                 reuse.pANN = F,
                                 sct = F)

  colnames(seurat_obj@meta.data)[grep("DF.classifications",colnames(seurat_obj@meta.data))] <- "DF.classfication1"
  seurat_obj <- doubletFinder(seurat_obj, PCs = 1:30, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = "DF.classfication1", sct = FALSE)
  colnames(seurat_obj@meta.data)[grep("DF.classifications_",colnames(seurat_obj@meta.data))] <- "DF.classfication2"
  return(seurat_obj)
}
