#==== library load in #####
pacman::p_load(ggplot2,patchwork,glue,dplyr)
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(EnsDb.Hsapiens.v86)
library(Matrix)
workdir <- "~"
setwd(workdir)
#==== Data load in #####
counts <- readMM("SNU601/peak/atac_matrix.mtx.gz")
CBs <- read.table("SNU601/peak/barcodes.tsv.gz")[,1]
peaks <- read.table("SNU601/peak/features.tsv.gz") |>
  mutate(peakid = paste0(V1,":",V2,"-",V3))
rownames(counts) <- peaks$peakid
colnames(counts) <- CBs
CNVclone <- readRDS("SNU601/clonality/SNU601_cloneinfo.rds")
counts.sub <- counts[,intersect(rownames(CNVclone),colnames(counts))]
CNVclone <- CNVclone[colnames(counts.sub),]
lowerDetect <-  0.1 * ncol(counts.sub)
set.seed(114)
#==== Create signac+seurat object #####
chrom_assay <- CreateChromatinAssay(
  counts = counts.sub,
  sep = c(":", "-"),
  genome = 'hg38',
  min.cells =lowerDetect,
  min.features = 0)
SNU601 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = CNVclone)
#==== Processing + clustering #####
SNU601 <- FindTopFeatures(SNU601, min.cells = 100)
SNU601 <- RunTFIDF(SNU601)
SNU601 <- RunSVD(SNU601, n = 100)
SNU601 <- RunUMAP(
  SNU601,
  reduction = "lsi",
  dims = 2:50,
  reduction.name = "UMAP")
SNU601 <- FindNeighbors(SNU601, dims = 2:50, reduction = "lsi")
SNU601 <- FindClusters(SNU601, resolution = 0.8, algorithm = 3)
trajInf<- as.cell_data_set(SNU601)
trajInf <- cluster_cells(cds = trajInf, reduction_method = "UMAP")
trajInf <- learn_graph(trajInf, use_partition = TRUE)
# rootclone <- FetchData(SNU601,c("UMAP_1","UMAP_2","sclone"))
# rootclone <-  rootclone |> dplyr::filter(UMAP_1< -1.8 & UMAP_2< -2 & sclone=="C1") 
# rootCB <- rownames(rootclone)
# saveRDS(rootCB,"SNU601/clonality/rootcb.rds")
#==== Pseudotime + trajectory analysis #####
rootCB <- readRDS("SNU601/clonality/rootcb.rds")
trajInf <- order_cells(trajInf, reduction_method = "UMAP", root_cells = rootCB)
SNU601 <- AddMetaData(
  object = SNU601,
  metadata = trajInf@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime")
saveRDS(SNU601,"SNU601/peak/SNU601Peak_Latest.rds")