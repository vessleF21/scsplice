# =============================================================================
# Doublet detection and removal using DoubletFinder
# Identifies and removes doublets (two cells captured in one droplet) from
# 10x Chromium scRNA-seq data. Uses a k-nearest neighbor approach to compute
# the proportion of artificial nearest neighbors (pANN) for each cell.
#
# Input:
#   --GSEname    : sample/dataset name
#   --Batchname  : batch identifier
#   --doublerate : expected doublet rate (e.g. 0.05 for 5%)
#
# Official DoubletFinder: https://github.com/chris-mcginnis-ucsf/DoubletFinder
# =============================================================================
library(optparse)

option_list = list(
  make_option("--GSEname",    type='character', help="Sample name"),
  make_option("--Batchname",  type='character', help="Batch name"),
  make_option("--doublerate", type='numeric',   help="Expected doublet rate")
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

# =============================================================================
# paramSweep: sweep pN and pK parameter space to find optimal pK
# pN: proportion of artificial doublets to simulate (5%-30%)
# pK: neighborhood size for pANN computation
# For datasets >10,000 cells, downsamples to 10,000 for efficiency
# =============================================================================
paramSweep <- function(seu, PCs=1:10, sct=FALSE, num.cores=1) {
  require(Seurat); require(fields); require(parallel)
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05,0.3,by=0.05)

  # Remove pK values that would produce fewer than 1 artificial doublet
  min.cells <- round(nrow(seu@meta.data)/(1-0.05) - nrow(seu@meta.data))
  pK <- pK[which(round(pK*min.cells) >= 1)]

  orig.commands <- seu@commands

  # Downsample to 10,000 cells if dataset is large
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000)]
  } else {
    real.cells <- rownames(seu@meta.data)
  }
  data <- seu@assays$RNA@counts[, real.cells]
  n.real.cells <- ncol(data)

  # Run parameter sweep in parallel or serial
  if (num.cores > 1) {
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)), FUN=parallel_paramSweep,
                        n.real.cells, real.cells, pK, pN, data,
                        orig.commands, PCs, sct, mc.cores=num.cores)
    stopCluster(cl)
  } else {
    output2 <- lapply(as.list(1:length(pN)), FUN=parallel_paramSweep,
                      n.real.cells, real.cells, pK, pN, data,
                      orig.commands, PCs, sct)
  }

  # Flatten nested list and assign pN_pK combination names
  sweep.res.list <- list()
  list.ind <- 0
  for (i in 1:length(output2))
    for (j in 1:length(output2[[i]])) {
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  name.vec <- NULL
  for (j in 1:length(pN))
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep="_"))
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}

# =============================================================================
# find.pK: identify optimal pK using bimodality coefficient (BCmvn)
# BCmvn = MeanBC / VarBC; higher value indicates better separation of
# real cells and artificial doublets in PCA space
# =============================================================================
find.pK <- function(sweep.stats) {
  '%ni%' <- Negate('%in%')
  if ("AUC" %ni% colnames(sweep.stats)) {
    bc.mvn <- as.data.frame(matrix(0L, nrow=length(unique(sweep.stats$pK)), ncol=5))
    colnames(bc.mvn) <- c("ParamID","pK","MeanBC","VarBC","BCmetric")
    bc.mvn$pK     <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)

    # Compute BCmvn for each pK across all pN values
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanBC[x]   <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x]    <- sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"]) /
                             sd(sweep.stats[ind, "BCreal"])^2
    }
    return(bc.mvn)
  }
}

# =============================================================================
# parallel_paramSweep: core function called per pN value
# 1. Simulates artificial doublets by averaging pairs of real cell profiles
# 2. Merges artificial doublets with real cells into a combined Seurat object
# 3. Runs standard preprocessing (normalize -> HVG -> scale -> PCA)
# 4. Computes pANN for each real cell across all pK values
#    pANN = proportion of k nearest neighbors that are artificial doublets
# =============================================================================
parallel_paramSweep <- function(n, n.real.cells, real.cells, pK, pN, data,
                                orig.commands, PCs, sct) {
  sweep.res.list <- list(); list.ind <- 0

  # Simulate artificial doublets by averaging random cell pairs
  n_doublets  <- round(n.real.cells/(1-pN[n]) - n.real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace=TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace=TRUE)
  doublets    <- (data[,real.cells1] + data[,real.cells2]) / 2
  colnames(doublets) <- paste("X", 1:n_doublets, sep="")
  data_wdoublets <- cbind(data, doublets)

  # Preprocess combined object (real + artificial doublets)
  if (!sct) {
    seu_wdoublets <- CreateSeuratObject(counts=data_wdoublets) %>%
      NormalizeData(normalization.method=orig.commands$NormalizeData.RNA@params$normalization.method,
                    scale.factor=orig.commands$NormalizeData.RNA@params$scale.factor) %>%
      FindVariableFeatures(selection.method=orig.commands$FindVariableFeatures.RNA$selection.method,
                           nfeatures=orig.commands$FindVariableFeatures.RNA$nfeatures) %>%
      ScaleData() %>%
      RunPCA(npcs=length(PCs), verbose=FALSE)
  } else {
    require(sctransform)
    seu_wdoublets <- CreateSeuratObject(counts=data_wdoublets) %>%
      SCTransform() %>%
      RunPCA(npcs=length(PCs))
  }

  # Compute PC distance matrix between all cells and real cells only
  nCells    <- nrow(seu_wdoublets@meta.data)
  pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[,PCs]
  rm(seu_wdoublets); gc()
  dist.mat  <- fields::rdist(pca.coord)[,1:n.real.cells]

  # Pre-order distances and trim to max pK neighborhood size
  for (i in 1:n.real.cells) dist.mat[,i] <- order(dist.mat[,i])
  dist.mat <- dist.mat[1:(round(nCells*max(pK))+5), ]

  # Compute pANN for each pK: fraction of k neighbors that are artificial doublets
  for (k in 1:length(pK)) {
    pk.temp <- round(nCells * pK[k])
    pANN <- data.frame(pANN=numeric(n.real.cells), row.names=real.cells)
    list.ind <- list.ind + 1
    for (i in 1:n.real.cells) {
      neighbors    <- dist.mat[2:(pk.temp+1), i]
      pANN$pANN[i] <- length(which(neighbors > n.real.cells)) / pk.temp
    }
    sweep.res.list[[list.ind]] <- pANN
  }
  return(sweep.res.list)
}

# =============================================================================
# removeDoublet: main wrapper function
# 1. Standard Seurat preprocessing (normalize, HVG, scale, PCA, UMAP, cluster)
# 2. Parameter sweep to find optimal pK via BCmvn
# 3. First DoubletFinder pass: classify doublets using expected doublet rate
# 4. Second DoubletFinder pass: adjust for homotypic doublets
#    (doublets from the same cell type, which are harder to detect)
#
# Output: Seurat object with two metadata columns:
#   DF.classfication1 : initial doublet classification
#   DF.classfication2 : adjusted classification accounting for homotypic doublets
# =============================================================================
removeDoublet <- function(seurat_obj) {

  # Standard preprocessing
  seurat_obj <- seurat_obj %>%
    NormalizeData(normalization.method="LogNormalize", scale.factor=10000) %>%
    FindVariableFeatures(selection.method="vst", nfeatures=3000) %>%
    ScaleData() %>%
    RunPCA(features=VariableFeatures(.), npcs=30) %>%
    RunUMAP(dims=1:30) %>%
    FindNeighbors() %>%
    FindClusters()

  # Find optimal pK via BCmvn metric
  seurat_test  <- paramSweep(seurat_obj, PCs=1:30, sct=FALSE)
  seurat_stats <- summarizeSweep(seurat_test, GT=FALSE)
  bcmvn        <- find.pK(seurat_stats)
  pK_bcmvn     <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  # Compute expected number of doublets from user-supplied doublet rate
  DoubletRate      <- as.numeric(opt$doublerate)
  homotypic.prop   <- modelHomotypic(seurat_obj$seurat_clusters)
  nExp_poi         <- round(DoubletRate * ncol(seurat_obj))
  nExp_poi.adj     <- round(nExp_poi * (1 - homotypic.prop))  # homotypic-adjusted count

  # First pass: classify using unadjusted expected doublet count
  seurat_obj <- doubletFinder(seurat_obj, pK=pK_bcmvn, PCs=1:30,
                               nExp=nExp_poi, reuse.pANN=FALSE, sct=FALSE)
  colnames(seurat_obj@meta.data)[grep("DF.classifications",
                                       colnames(seurat_obj@meta.data))] <- "DF.classfication1"

  # Second pass: reclassify using homotypic-adjusted doublet count
  seurat_obj <- doubletFinder(seurat_obj, PCs=1:30, pN=0.25, pK=pK_bcmvn,
                               nExp=nExp_poi.adj, reuse.pANN="DF.classfication1", sct=FALSE)
  colnames(seurat_obj@meta.data)[grep("DF.classifications_",
                                       colnames(seurat_obj@meta.data))] <- "DF.classfication2"
  return(seurat_obj)
}





