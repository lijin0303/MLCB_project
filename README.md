# Improved variant calling and phylogeny inference from scATAC-seq

Scripts, results and documentations for MLCB group project 2023

Megan Le, Ruitong Li, Daniel Schaffer

Single-cell, genome-wide profiling of accessible chromatin (scATAC-seq) grants unprecedented scale and resolution to dissect phenotypic heterogeneity in the regulatory landscape. However, understanding changes in gene regulation along clonal evolution and linking them to any individual/disease phenotype remains a major open problem.
There have been emerging efforts in calling SNVs from sc-omics sequencing data (mainly scRNA-seq) in cancer samples. But they either suffer badly from poor accuracy or require specific sample composition such as diverse normal cell types in combination with tumor cells) (Muyas et al. 2023).

Compared to scRNA-seq, there are several advantages to using scATAC-seq data. First, scATAC-seq can be more reliably performed on FFPE/aged samples. Additionally, ATACseq covers select non coding regions of the genome that are never sequenced in RNAseq, and the coding regions of the genome are (presumably) under more constraint. Finally, RNA editing events may introduce spurious mutations in RNA-seq data. Therefore, mutation calling from profiles of accessible regions instead of mRNA abundance is preferred and will theoretically serve as a better lineage marker. 
Mutation calling from scATAC-seq is still a field left largely unexplored. (Dou et al. 2023) cleverly used phasing of SNPs to improve the accuracy of nuclear single nucleotide variants (SNVs), which are found to be consistent with lineage information (eg: differentiation during hematopoiesis). This method is called Monopgen.
