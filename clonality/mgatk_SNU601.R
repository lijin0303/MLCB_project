setwd("~/SNU601/")
source("~/InHouse/Basics.R")
pacman::p_load(ComplexHeatmap,InteractiveComplexHeatmap,
               RColorBrewer,wesanderson,
               usedist,circlize)
# mgatk <- read.csv("SNU601_mgatk_het.csv",check.names = F,row.names  = 1)
mgatk <- data.table::fread("mgatk.cell_heteroplasmic_df.tsv.gz") |> 
  column_to_rownames("V1")
cnvClone <- readRDS("umap_clone.rds")
cell2focus <- intersect(rownames(cnvClone),rownames(mgatk))
cell2focus <- cnvClone[cell2focus,] |> 
  rownames_to_column("CB") |> 
  arrange(Clone) |> 
  pull(CB)
# mgatkrds <- readRDS("mgatk.rds")
# SummarizedExperiment::assays(mgatkrds)[["A_counts_fw"]]
#==== Heatmap: cell by mut =====
col_fun = colorRamp2(c(0,0.5,1),c("beige","orange","brown"))
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
clone2focus <- sort(unique(cnvClone[cell2focus,"Clone"]))
clonecol = brewer.pal(n = length(clone2focus), name = "Dark2")
names(clonecol) = clone2focus
rowannot = rowAnnotation(
  col = list(clone = clonecol),
  clone = cnvClone[cell2focus,"Clone"],
  gp = gpar(col = "black",lwd = 0.01))
cm_cluster <- Heatmap(matrix = mgatk[cell2focus,],name="mtSNV",
                      rect_gp = gpar(col = "black", lwd = 0.1),
                      na_col = "black",
                      col = col_fun,
                      ## column setting ##
                      show_column_names = F,
                      column_names_rot = 45, 
                      column_names_side = "top",
                      
                      cluster_columns = T,
                      clustering_distance_columns = jaccard,
                      show_column_dend = T,
                      column_dend_side = "top",
                      column_dend_height = unit(2, "cm"),
                      
                      ## row setting ##
                      show_row_names = F,
                      cluster_rows = F,
                      show_row_dend = T,
                      clustering_distance_rows = jaccard,
                      left_annotation=rowannot,
                      show_heatmap_legend=T)

ht_shiny(cm_cluster, width1=500,height1=600)
