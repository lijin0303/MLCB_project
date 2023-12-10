source("~/InHouse/dar_utils.R")
require(RColorBrewer)
jaccard <- function(x1, x2){
  ax1 = x1
  ax2 = x2
  ax1[is.na(ax1)] = 0
  ax2[is.na(ax2)] = 0
  na.id <- is.na(x1) | is.na(x2)
  alt <- ax1!=0 | ax2!=0
  denom = sum((!na.id) & alt)
  if(denom>0){
    jd <- 1- (sum(ax1!=0 & ax2!=0)/denom)
  }else{
    jd <- 0.95
  }
  return(jd)
}
#=== Figure 1: somatic mut filtering ######
inputP <- "~/SNU601/clonality/scATAC"
cellname <- read.table(glue("{inputP}/chr1.cell.txt"))
CBmut_load <- map(1:22,\(chr){
  message("chr",chr)
  SNVmat <- readRDS(glue("{inputP}/chr{chr}.SNV_mat.RDS"))
  CB_SNV <- SNVmat[,-c(1:18)] |> t()
  CBmut <- CB_SNV[,apply(CB_SNV,2,\(x) sum(x!="0/0"))>0.05*nrow(cellname),drop=F]
  return(CBmut)
}) |> bind_cols()
CBmut_all <- as.matrix(CBmut_load)
rownames(CBmut_all) <- cellname$V1
rownames(CBmut_load) <- cellname$V1
CBmut_all <- CBmut_all[apply(CBmut_all,1,\(x) sum(x!="0/0"))>0.1*ncol(CBmut_all),,drop=F]
CBmut_bi <- CBmut_all
CBmut_bi[CBmut_bi=="0/0"] <- "ND"
CBmut_bi[CBmut_bi=="1/0"] <- "WT"
CBmut_bi[CBmut_bi=="1/1"] <- "ALT"
CBmut_bi[CBmut_bi=="0/1"] <- "ALT"
cellno <- nrow(CBmut_bi)
var_summary <- data.frame(af = apply(CBmut_bi,2,\(x) sum(x=="ALT"))/nrow(CBmut_bi))
require(lessR)
Histogram(af,data=var_summary,
          xlab = "scATAC-seq: Alternative Allele Frequency",
          ylab = "Number of Variants")
pl <- 0.1
ph <- 0.9
CBmut_bi <- CBmut_bi[,apply(CBmut_bi,2,\(x) sum(x=="ALT"))>pl*cellno,drop=F]
CBmut_bi <- CBmut_bi[,apply(CBmut_bi,2,\(x) sum(x=="ALT"))<ph*cellno,drop=F]
CB_clone <- readRDS("~/SNU601/clonality/umap_clone.rds") |>
  mutate(Clone=case_when(Clone==3~4,
                         Clone==4~3,
                         TRUE~Clone)) |>  # the 1-6 defined in umap-clone alleloscope is wrongly ordered by dendrogram
  mutate(clone=paste0("clone_",Clone)) |> 
  rownames_to_column("CB") |> 
  mutate(sclone = paste0("C",Clone))
CB_cloneSort <- CB_clone%>%arrange(clone) |> column_to_rownames("CB")
saveRDS(CB_clone%>%arrange(sclone) |> column_to_rownames("CB"),"~/SNU601/clonality/SNU601_cloneinfo.rds")

CB_cloneSort <- CB_cloneSort[
  intersect(rownames(CB_cloneSort),rownames(CBmut_bi)),,drop=F]
clonemap <-CB_cloneSort$clone
#clonecol = brewer.pal(n = 6, name = "Dark2")
clonecol = c("#CD3217", "#30A621", "#E97A0F", "#568ABC", "#843BB6", "#B85F24")
names(clonecol) = paste0("clone_",1:6)
rowannot = rowAnnotation(
  col = list(clone = clonecol),
  clone = clonemap,
  gp = gpar(col = "black",lwd = 0.01))
numericMat <- CBmut_bi
numericMat[numericMat=="ND"] <- NA
numericMat[numericMat=="ALT"] <- 1
numericMat[numericMat=="WT"] <- 0
mat2col = c(`2`="#63bff0",`1`="#9D0208",`0`="gray")
numericMat <- matrix(as.numeric(numericMat), ncol = ncol(numericMat))  
colnames(numericMat) <- colnames(CBmut_bi)
rownames(numericMat) <- rownames(CBmut_bi)
mat2plt <- numericMat[rownames(CB_cloneSort),]
saveRDS(list("ATAC_call"=CBmut_load,"mat2plt"=mat2plt,
             "rowannot"=rowannot,"matcol"=mat2col,
             "AFsummary"=var_summary),"~/SNU601/clonality/scATAC_mutcall.rds")
var2zoom <- apply(mat2plt,2,\(x) table(x))
refalt <- apply(var2zoom,2,\(c) sum(c>40)==2) # ref & var detected in >20 cells
sumcov <- colSums(var2zoom)>200
mut_topcov <- mat2plt[,refalt&sumcov]
saveRDS(mut_topcov,"~/SNU601/clonality/scATAC_highcovmut.rds")
#=== Figure 2: Heatmap plotting ######
somut <- readRDS("~/SNU601/clonality/scATAC_mutcall.rds")
scmutCov <- data.table::fread("~/SNU601/cxm_cov_mat.txt")
scmutCov <-  as.matrix(scmutCov)
rownames(scmutCov) <- rownames(somut$ATAC_call)
colnames(scmutCov) <- colnames(somut$ATAC_call)
cellcled = Heatmap(matrix = somut$mat2plt,name=" ",
                   rect_gp = gpar(col = "white", lwd = 0.05),
                   na_col = "white",
                   col = somut$matcol,
                   ## column setting ##
                   show_column_names = F,
                   cluster_columns = T,
                   clustering_distance_columns = jaccard,
                   clustering_distance_rows= jaccard,
                   ## row setting ##
                   show_row_names = F,
                   cluster_rows = T,
                   show_row_dend = T,
                   show_heatmap_legend=T,
                   left_annotation=somut$rowannot)
CNVcled = Heatmap(matrix = somut$mat2plt,name=" ",
                   rect_gp = gpar(col = "white", lwd = 0.05),
                   na_col = "white",
                   col = somut$matcol,
                   ## column setting ##
                   show_column_names = F,
                   cluster_columns = T,
                   clustering_distance_columns = jaccard,
                   clustering_distance_rows= jaccard,
                   ## row setting ##
                   show_row_names = F,
                   cluster_rows = F,
                   show_row_dend = T,
                   show_heatmap_legend=T,
                   left_annotation=somut$rowannot)

covmat <- scmutCov[rownames(somut$mat2plt),column_order(CNVcled)]
quantile(covmat,c(0.9,1,0.01))
# covmat <- covmat[1:20,1:30]
# library(circlize)
# pl = list()
# for(palette in sort(hcl.pals("sequential"))) {
#   palette = "Rocket"
#   col_fun = colorRamp2(seq(0,10,1), hcl_palette = palette)
#   # ht = Heatmap(mat, name = "mat", col = col_fun, 
#   #              show_row_dend = FALSE, show_column_dend = FALSE,
#   #              column_title = paste0("palette = '", palette, "'"),
#   #              heatmap_legend_param = list(legend_height = unit(6, "cm")))
#   covmat <- scmutCov[rownames(somut$mat2plt),column_order(CNVcled)]
#   covmat[covmat==0] <- NA
#   coverage = Heatmap(matrix = covmat,name=" ",
#                      rect_gp = gpar(col = "white", lwd = 0.05),
#                      na_col = "white",
#                      col = col_fun, 
#                      show_row_dend = FALSE, show_column_dend = FALSE,
#                      column_title = paste0("palette = '", palette, "'"),
#                      ## column setting ##
#                      show_column_names = F,
#                      cluster_columns = F,
#                      ## row setting ##
#                      show_row_names = F,
#                      cluster_rows = F,
#                      heatmap_legend_param = list(legend_height = unit(6, "cm")))
#   coverage
#   pl[[palette]] = grid.grabExpr(draw(coverage))
# }
# library(cowplot)
# plot_grid(plotlist = pl, ncol = 5)


saveRDS(list("rowcl"=cellcled,"cloneOrd"=CNVcled),"~/SNU601/clonality/heatmap_mutmat.rds")

png("~/SNU601/CellClustered_cellmutmat.png",
    width=5.5,height=8.5,units = "in",res=1200)
draw(htMut$rowcl)
dev.off()

png("~/SNU601/CloneOrdered_cellmutmat.png",
    width=5.5,height=8.5,units = "in",res=1200)
draw(htMut$cloneOrd)
dev.off()
# pdf("~/SNU601/CloneOrdered_cellmutmat_coverage.pdf",width = 5,height=5)
# draw(coverage)
# dev.off()
# 
# png("~/SNU601/CellClustered_cellmutmat.png",res=500)
# draw(cellcled)
# dev.off()

ht_shiny(cellcled)
ht_shiny(CNVcled)
#==== Reference #####
setwd("~/SNU601/")
source("~/InHouse/Basics.R")
chr = 3
cellname <- read.table(glue("~/SNU601/chr{chr}.cell.txt"))
SNVmat <- readRDS(glue("~/SNU601/chr{chr}.SNV_mat.RDS"))
CB_SNV <- SNVmat[,-c(1:18)] |> t()
rownames(CB_SNV) <- cellname$V1
min_cov = round(nrow(cellname)*0.1,0)
CBmut <- CB_SNV[,apply(CB_SNV,2,\(x) sum(x!="0/0"))>min_cov]



biCBmut <- CBmut
biCBmut[biCBmut!="1/0"] <- ""
cnvClone[cell2focus,"Clone"]
apply(CBmut[cell2focus,],2,\(x)table(cnvClone[cell2focus,"Clone"],x))




mat2plt <- CBmut
cnvClone <- readRDS("clonality/umap_clone.rds")
cell2focus <- intersect(rownames(cnvClone),rownames(mat2plt))
cell2focus <- cnvClone[cell2focus,] |> 
  rownames_to_column("CB") |> 
  arrange(Clone) |> 
  pull(CB)
library(RColorBrewer)
library(ComplexHeatmap)
clone2focus <- sort(unique(cnvClone[cell2focus,"Clone"]))
clonecol = brewer.pal(n = length(clone2focus), name = "Dark2")
names(clonecol) = clone2focus
rowannot = rowAnnotation(
  col = list(clone = clonecol),
  clone = cnvClone[cell2focus,"Clone"],
  gp = gpar(col = "black",lwd = 0.01))


mat2col = c("white","#63bff0","#9D0208","gray")
names(mat2col) = c("0/0","0/1","1/1","1/0")

jaccard <- function(x1, x2){
  ax1 = x1
  ax2 = x2
  ax1[is.na(ax1)] = 0
  ax2[is.na(ax2)] = 0
  na.id <- is.na(x1) | is.na(x2)
  alt <- ax1!=0 | ax2!=0
  denom = sum((!na.id) & alt)
  if(denom>0){
    jd <- 1- (sum(ax1!=0 & ax2!=0)/denom)
  }else{
    jd <- 0.95
  }
  return(jd)
}

numericMat <- mat2plt
numericMat[numericMat=="0/0"] <- NA
numericMat[numericMat=="0/1"] <- 2
numericMat[numericMat=="1/0"] <- 0
numericMat[numericMat=="1/1"] <- 1
mat2col = c(`2`="#63bff0",`1`="#9D0208",`0`="gray")
numericMat <- matrix( as.numeric(numericMat), ncol = ncol(numericMat))  
colnames(numericMat) <- colnames(mat2plt)
rownames(numericMat) <- rownames(mat2plt)
# mat2plt <- numericMat[cell2focus,1:30]
ht2check = Heatmap(matrix = numericMat[cell2focus,],name=" ",
                   rect_gp = gpar(col = "white", lwd = 0.05),
                   na_col = "white",
                   col = mat2col,
                   ## column setting ##
                   show_column_names = F,
                   cluster_columns = T,
                   clustering_distance_columns = jaccard,
                   clustering_distance_rows= jaccard,
                   ## row setting ##
                   show_row_names = F,
                   cluster_rows = T,
                   show_row_dend = T,
                   show_heatmap_legend=T,
                   left_annotation=rowannot)
library(InteractiveComplexHeatmap)
ht_shiny(ht2check)
