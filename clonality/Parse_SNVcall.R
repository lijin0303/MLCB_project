chr = 1
cellname <- read.table(glue("~/SNU601/chr{chr}.cell.txt"))
SNVmat <- readRDS(glue("~/SNU601/chr{chr}.SNV_mat.RDS"))
CB_SNV <- SNVmat[,-c(1:18)] |> t()
rownames(CB_SNV) <- cellname$V1
min_cov = round(nrow(cellname)*0.1,0)
CBmut <- CB_SNV[,apply(CB_SNV,2,\(x) sum(x!="0/0"))>min_cov]

mat2plt <- CBmut
cnvClone <- readRDS("umap_clone.rds")
cell2focus <- intersect(rownames(cnvClone),rownames(mat2plt))
cell2focus <- cnvClone[cell2focus,] |> 
  rownames_to_column("CB") |> 
  arrange(Clone) |> 
  pull(CB)
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
    jd <- 0.5
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
mat2plt <- numericMat[cell2focus,1:30]
ht2check = Heatmap(matrix = numericMat[cell2focus,1:30],name=" ",
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
                   cluster_rows = F,
                   show_row_dend = T,
                   show_heatmap_legend=T,
                   left_annotation=rowannot)
ht_shiny(ht2check)
