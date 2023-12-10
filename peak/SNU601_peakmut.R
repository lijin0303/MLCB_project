pacman::p_load(ggplot2,patchwork,glue,dplyr,tidyr,magrittr,ggpubr)
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(EnsDb.Hsapiens.v86)
library(Matrix)
workdir <- "~"
setwd(workdir)
SNU601_peak <- readRDS("SNU601/peak/SNU601Peak_Latest.rds")
SNU601_mut <- readRDS("SNU601/peak/SNU601Mut_Latest.rds")
SNU601_peak <- AddMetaData(
  object = SNU601_peak,
  metadata = SNU601_mut$pseudotime,
  col.name = "pseudotimeMut")
png("~/SNU601/umap_PeakCluster.png",
    width=5.5,height=5,units = "in",res=1200)
colors_ <- ggsci::pal_npg()(10)
DimPlot(SNU601_peak, 
        label = TRUE,label.color = "white",
        label.box = T,label.size=rel(3))+
  scale_color_manual(values = colors_ )+
  scale_fill_manual(values = colors_ )+
  labs(x="UMAP 1",y="UMAP 2")+
  theme(aspect.ratio = 1,
        axis.title = element_text(size = rel(0.6)))+
  ggtitle("Colored by scATAC \n clusters with LSI")
dev.off()


png("~/SNU601/umap_mutclone.png",
    width=5.5,height=5,units = "in",res=1200)
clonecol <- c("#CD3217", "#30A621", "#E97A0F", "#568ABC", "#843BB6", "#B85F24")
names(clonecol) <- paste0("C",1:6)
DimPlot(SNU601_peak, group.by="sclone",
        label = TRUE,label.color = "white",
        label.box = T,label.size=rel(3))+
  scale_color_manual(values = clonecol)+
  scale_fill_manual(values = clonecol)+
  labs(x="UMAP 1",y="UMAP 2")+
  theme(aspect.ratio = 1,axis.title = element_text(size = rel(0.6)))+
  ggtitle("Colored by CNV Clones")
dev.off()

png("~/SNU601/umap_PeakPseudo.png",
    width=5.5,height=5,units = "in",res=1200)
FeaturePlot(SNU601_peak,features = 'pseudotime', pt.size = 0.1) +
  labs(x="UMAP 1",y="UMAP 2")+
  theme(aspect.ratio = 1,axis.title = element_text(size = rel(0.6)))+ 
  scale_color_viridis_c()
dev.off()

g1 <- FeaturePlot(SNU601_peak,features = 'pseudotime', pt.size = 0.1) +
  labs(x="UMAP 1",y="UMAP 2")+
  theme(aspect.ratio = 1,axis.title = element_text(size = rel(0.6)))+ 
  scale_color_viridis_c()
g2 <- FeaturePlot(SNU601_peak,features = 'pseudotimeMut', pt.size = 0.1) +
  labs(x="UMAP 1",y="UMAP 2",title="Somatic Mutation-based")+
  theme(aspect.ratio = 1,axis.title = element_text(size = rel(0.6)))+ 
  scale_color_viridis_c()
pseudoTime <- ggpubr::ggarrange(g1,g2,nrow=1)
png("~/SNU601/umap_PeakPseudo.png",
    width=9,height=5,units = "in",res=1200)
pseudoTime
dev.off()
#==== Peak signal in Ref vs Alt ====
mut_topcov <-  readRDS("~/SNU601/clonality/scATAC_highcovmut.rds")
SNU601 <- readRDS("~/SNU601/peak/SNU601_clustered.rds")
mutgr <- data.frame(varName= colnames(mut_topcov)) |> 
  separate(varName,into = c("chr","start","ref","alt"),sep=":") |> 
  mutate(end=start) |> 
  makeGRangesFromDataFrame(keep.extra.columns = T)
mutcontext <- ArchR::extendGR(gr = mutgr, upstream = 50000, downstream = 50000)
mutregions <- paste0(seqnames(mutcontext),"-",ranges(mutcontext))
peak_nearby <- findOverlaps(mutcontext,granges(SNU601@assays$peaks))
peakDist <- function(i){
  nonzero_peakSignal <- mut_topcov[,i,drop=F] |> 
    set_colnames(c("mut")) |> 
    cbind(SNU601@assays$peaks$data[subjectHits(peak_nearby[queryHits(peak_nearby)==i,]),
                                   rownames(mut_topcov)] |> 
            as.matrix() |> t()) |> 
    as.data.frame() |> 
    mutate(mut = case_when(mut==1~"Alt",
                           mut==0~"Ref",
                           TRUE~"ND")) |>
    gather(peak,cnt,-mut) |> 
    filter(cnt>0) 
  p1 <- ggplot(nonzero_peakSignal, aes(x=cnt, fill=mut)) +
    geom_histogram(color="#e9ecef", alpha=0.8,binwidth=0.1) +
    scale_fill_manual(values=c("#69b3a2", "#404080","orange")) +
    theme_bw() +
    labs(fill="",x="TF-IDF normalized peak signal",y="Frequency")+
    guides(fill="none")
  my_comparisons <- list( c("Alt", "Ref"))
  p2 <- ggplot(nonzero_peakSignal, aes(x = mut, y = cnt, fill=mut))+
    scale_fill_manual(values=c("#69b3a2", "#404080","orange")) + 
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1, fill="white")+
    theme_minimal() +
    labs(fill="",x="",y="Normalized Peak Signal")+ 
    stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = 2.5)     # Add global p-value
  pcomb <- ggarrange(p1,p2,nrow = 1)
  pcomb <- annotate_figure(pcomb, top = text_grob(glue("Mutation: {colnames(mut_topcov)[i]}"), 
                                                  face = "bold", size = rel(10)))
  return(pcomb)
}
peakDistOuts <- purrr::map(1:ncol(mut_topcov),peakDist)
DA_result <- purrr::map(1:ncol(mut_topcov),\(x) {
  singlemut = case_when(mut_topcov[,i]==1~"Alt",mut_topcov[,i]==0~"Ref",TRUE~"ND") 
  names(singlemut) <- rownames(mut_topcov)
  SNU601 <- AddMetaData(
    object = SNU601,
    metadata = singlemut,
    col.name = "mutstatus")
  da_peaks <- FindMarkers(
    object = SNU601,
    group.by="mutstatus",
    ident.1 = "Alt",
    ident.2 = "Ref",
    test.use = 'wilcox')
  print(nrow(da_peaks))
  da_peaks%<>%slice(1:5)
  return(da_peaks)
})
for(i in 1:34){
  print(i)
  png(glue("~/SNU601/peakDist/{colnames(mut_topcov)[i]}_nearbyPeak.png"),
      width=3500,height=1700,res=300)
  print(peakDistOuts[[i]])
  dev.off() 
}

