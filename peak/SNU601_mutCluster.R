pacman::p_load(ggplot2,patchwork,glue,dplyr,purrr)
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(EnsDb.Hsapiens.v86)
library(Matrix)
workdir <- "~"
setwd(workdir)
somut <- readRDS("~/SNU601/clonality/scATAC_mutcall.rds")
counts <- somut$ATAC_call
counts[counts=="0/0"] <- "0"
counts[counts=="1/0"] <- "0"
counts[counts=="1/1"] <- "1"
counts[counts=="0/1"] <- "1"
numericMat <- as.matrix(counts)
numericMat <- matrix(as.numeric(numericMat), ncol = ncol(numericMat))  
colnames(numericMat) <- colnames(somut$ATAC_call)
rownames(numericMat) <- rownames(somut$ATAC_call)
numericMat <- t(numericMat)
CNVclone <- readRDS("SNU601/clonality/SNU601_cloneinfo.rds")
counts.sub <- numericMat[,intersect(rownames(CNVclone),colnames(numericMat))]
CNVclone <- CNVclone[colnames(counts.sub),]
lowerDetect <-  0.01 * ncol(counts.sub)
rownames(counts.sub) <- map_chr(rownames(counts.sub) ,\(x){
  grstr <- strsplit(x,split=":")[[1]][c(1,2,2)]
  grstr <- paste0(grstr[1],":",grstr[2],"-",grstr[3])
  return(grstr)
})
set.seed(114)
mut_assay <- CreateChromatinAssay(
  counts = counts.sub,
  sep = c(":","-"),
  genome = 'hg38',
  min.cells =lowerDetect,
  min.features = 0)
SNU601 <- CreateSeuratObject(
  counts = mut_assay,
  assay = "mut",
  meta.data = CNVclone)
SNU601 <- FindTopFeatures(SNU601, min.cells = 5)
SNU601 <- RunTFIDF(SNU601)
SNU601 <- RunSVD(SNU601, n = 100)
SNU601 <- RunUMAP(
  SNU601,
  reduction = "lsi",
  dims = 1:50,
  reduction.name = "UMAP")
SNU601 <- FindNeighbors(SNU601, dims = 1:50, reduction = "lsi")
SNU601 <- FindClusters(SNU601, resolution = 0.85, algorithm = 3)
trajInf<- as.cell_data_set(SNU601)
trajInf <- cluster_cells(cds = trajInf, reduction_method = "UMAP")
trajInf <- learn_graph(trajInf, use_partition = TRUE)
# rootclone <- FetchData(SNU601,c("UMAP_1","UMAP_2","sclone"))
# rootclone <-  rootclone |> dplyr::filter(UMAP_1< -1.8 & UMAP_2< -2 & sclone=="C1") 
# rootCB <- rownames(rootclone)
# saveRDS(rootCB,"SNU601/clonality/rootcb.rds")
rootCB <- readRDS("SNU601/clonality/rootcb.rds")
trajInf <- order_cells(trajInf, reduction_method = "UMAP", root_cells = rootCB)
SNU601 <- AddMetaData(
  object = SNU601,
  metadata = trajInf@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime")
saveRDS(SNU601,"SNU601/peak/SNU601Mut_Latest.rds")