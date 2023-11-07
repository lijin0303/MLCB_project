library(Alleloscope)
setwd("~/MLCB/")
options(stringsAsFactors=F)
dir_path <- "SNU601/output";dir.create(dir_path)
inputF <- "SNU601/scATAC"
pacman::p_load(glue,dplyr)
data(centromere.GRCh38)
data(telomere.GRCh38)
size=read.table("sizes.cellranger-GRCh38-1.0.0.txt", 
                stringsAsFactors = F)
# SNP by cell matrices for ref and alt alleles
barcodes=read.table("SNU601/scATAC/barcodes.tsv", sep='\t', 
                    stringsAsFactors = F, header=F)
alt_all=readMM("SNU601/scATAC/alt_all.mtx")
ref_all=readMM("SNU601/scATAC/ref_all.mtx")
var_all=read.table("SNU601/scATAC/var_all.vcf", header = F, 
                   sep='\t', stringsAsFactors = F)
# bin by cell matrices for tumor and normal for segmentation
raw_counts=read.table('SNU601/scATAC/chr200k_fragments_sub.txt', 
                      sep='\t', header=T, row.names = 1,stringsAsFactors = F)
Obj=Createobj(alt_all =alt_all, ref_all = ref_all, var_all = var_all ,
              samplename='SNU601',
              genome_assembly="GRCh38",
              dir_path=dir_path, barcodes=barcodes, 
              size=size, assay='scATACseq')
Obj_filtered=Matrix_filter(Obj=Obj, cell_filter=5, SNP_filter=5, 
                           centro=centromere.GRCh38, 
                           telo=telomere.GRCh38) 
Obj_scDNA=readRDS("SNU601/scATAC/SNU601_dna.rds")

Obj_filtered$seg_table_filtered=Obj_scDNA$seg_table_filtered
Obj_filtered = Est_regions(Obj_filtered = Obj_filtered, 
                           max_nSNP = 30000, min_cell = 20, 
                           phases = Obj_scDNA$rds_list, plot_stat = T, 
                           cont = TRUE)
Obj_filtered$ref=Obj_scDNA$ref
Obj_filtered=Genotype_value(Obj_filtered = Obj_filtered, type='cellline', 
                            raw_counts=raw_counts, cov_adj =1,
                            ref_gtv = Obj_scDNA$genotype_values) 
Obj_filtered=Genotype(Obj_filtered = Obj_filtered, 
                               ref_gt = Obj_scDNA$genotypes,xmax=4)
clone.genotypes=readRDS("SNU601/scATAC/clone.genotypes.rds")
Obj_filtered=AssignClones_ref(Obj_filtered=Obj_filtered, clone.genotypes=clone.genotypes)
umap_peak=readRDS("SNU601/scATAC/peak_umap.rds")
Clone=Obj_filtered$cloneAssign$cloneAssign[match(rownames(umap_peak), names(Obj_filtered$cloneAssign$cloneAssign))]
umap_peak=cbind(umap_peak, Clone)