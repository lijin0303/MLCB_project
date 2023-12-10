#==== Data loadin ====
library(MutationalPatterns)
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(gridExtra)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
myplot_spectrum_region<- function (type_occurrences) {
  COLORS6 <- c(
    "#2EBAED", "#000000", "#DE1C14",
    "#D4D2D2", "#ADCC54", "#F0D0CE"
  )
  `C>T at CpG` <- `C>T other` <- type <- amount <- stdev <- tot_muts <- lower <- upper <- freq <- NULL
  total_indv <- sem <- error_95 <- NULL
  mode <- "relative_sample_feature"
  row_names <- rownames(type_occurrences)
  max_dots_in_name <- row_names %>% stringr::str_count("\\.") %>% 
    max()
  sample_names <- stringr::str_remove(row_names, "\\..*")
  feature <- stringr::str_remove(row_names, ".*\\.")
  feature <- factor(feature, levels = unique(feature))
  type_occurrences <- type_occurrences %>% dplyr::select(-`C>T at CpG`, 
                                                         -`C>T other`)
  tb_per_sample <- type_occurrences %>% as.data.frame() %>% 
    tibble::as_tibble() %>% dplyr::mutate(sample = sample_names, 
                                          feature = feature) %>% 
    tidyr::gather(key = "type", value = "amount", -sample, -feature)
  tb_per_sample <- tb_per_sample %>% dplyr::group_by(sample, 
                                                     feature) %>% dplyr::mutate(freq = amount/sum(amount)) %>% 
    dplyr::ungroup()
  y_lab <- "Relative contribution"
  tb_per_sample <- dplyr::mutate(tb_per_sample, freq = ifelse(is.nan(freq), 0, freq))
  by <- "all"
  tb_by <- tibble::tibble(sample = unique(tb_per_sample$sample), by = by)
  tb_per_sample <- tb_per_sample %>% dplyr::left_join(tb_by,by = "sample")
  tb <- tb_per_sample %>% dplyr::mutate(by = factor(by, levels = unique(by))) %>% 
    dplyr::group_by(by, feature, type) %>% 
    dplyr::summarise(stdev = sd(freq), freq = mean(freq),
                     amount = sum(amount), total_indv = dplyr::n(), 
                     .groups = "drop_last") %>% 
    dplyr::ungroup()
  spacing <- 0.5
  fig <- ggplot(tb, aes(x = feature, y = freq, fill = type, alpha = feature)) + 
    geom_bar(stat = "identity", position = "dodge", colour = "black", cex = 0.5) + 
    facet_grid(. ~ type) +
    scale_fill_manual(values = COLORS6) + 
    scale_alpha_discrete(range = c(1, 0.4)) + 
    labs(y = y_lab, x = "",fill="") + 
    theme_bw() +
    scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    theme(panel.spacing.x = unit(spacing, "lines"))+ 
    guides(alpha = "none")
  return(fig)
}
scATACmut <- readRDS("~/SNU601/clonality/scATAC_mutcall.rds")
mut5pcov <- colnames(scATACmut$ATAC_call)
options(stringsAsFactors = F)
pacman::p_load(optparse,glue,mgsub,dplyr,pbmcapply,pbapply,ggplot2,lessR)
pacman::p_load(rtracklayer,purrr)
import::from(tidyr,unite,separate,gather,spread)
import::from(magrittr,set_colnames,set_rownames,"%<>%")
import::from(tibble,column_to_rownames,rownames_to_column)
#==== Mutation signature ====
mutdf <- data.frame(mutcol=mut5pcov) |> 
  separate(mutcol,into=c("chr","start","REF","ALT"),sep=":") |>
  mutate(end=start)
mutgr <- makeGRangesFromDataFrame(mutdf,keep.extra.columns = T)
names(mutgr) <- gsub("chr","",paste0(mutdf$chr,":",mutdf$start,"_",mutdf$REF,"/",mutdf$ALT))
muts <- mutations_from_vcf(mutgr)
GenomeInfoDb::genome(mutgr) = 'hg38'
mut_context(mutgr,ref_genome)
mut_mat <- mut_matrix(vcf_list = mutgr, ref_genome = ref_genome)
sig96_mut <- plot_96_profile(mut_mat)
signatures = get_known_signatures()
fit_res <- fit_to_signatures(mut_mat, signatures)
knownsig <- plot_contribution(fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "absolute")
#==== Strand bias: transcription ====
genes_hg38 <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
mut_mat_s <- mut_matrix_stranded(mutgr, ref_genome, genes_hg38)
plot_192_profile(mut_mat_s)
strand_counts <- strand_occurrences(mut_mat_s)
strand_bias <- strand_bias_test(strand_counts)
ps1 <- plot_strand(strand_counts, mode = "relative")
ps2 <- plot_strand_bias(strand_bias, sig_type = "p")
strand_bias_notstrict <- strand_bias_test(strand_counts,
                                          p_cutoffs = c(0.5, 0.1, 0.05),
                                          fdr_cutoffs = 0.5)
#==== Genomic distribution ====
chromosomes <- seqnames(get(ref_genome))[1:22]
g <- getBSgenome(ref_genome, masked=FALSE)
seqlengths(mutgr) <- seqlengths(g)[1:22]
# Make a rainfall plot
CTCF_g <- readRDS(system.file("states/CTCF_g_data.rds",
                              package = "MutationalPatterns"))
promoter_g <- readRDS(system.file("states/promoter_g_data.rds",
                                  package = "MutationalPatterns"))
flanking_g <- readRDS(system.file("states/promoter_flanking_g_data.rds",
                                  package = "MutationalPatterns"))
regions <- GRangesList(promoter_g, flanking_g, CTCF_g)
ch = import.chain("~/SNU601/hg19ToHg38.over.chain")
names(regions) <- c("Promoter", "Promoter flanking", "CTCF")
seqlevelsStyle(regions) <- "UCSC"
regions <- map(regions,\(x) unlist(liftOver(x, ch)))
listToGRangesList = function(lst) {
  if(!is(lst, "GRangesList")) {
    if ("list" %in% class(lst)) {
      #strip elementMetadata
      lst = lapply(lst, function(x) { values(x) <- NULL; return(x) } )
      lst = GRangesList(lst)
    } else {
      message(cleanws("Converting GRanges to GRangesList."))
      lst = GRangesList(list(lst))
    }
  }
  return(lst)
}
regions <- listToGRangesList(regions)
grl_region <- split_muts_region(mutgr, regions)
type_occurrences_region <- mut_type_occurrences(grl_region, ref_genome)
#==== All visualization #####
g1 <- plot_192_profile(mut_mat_s)+
  theme(legend.position = "top")
g2 <- plot_strand_bias(strand_bias_notstrict, sig_type = "p")
strandbias <- ggpubr::ggarrange(g1,g2,ncol=1)
png("~/SNU601/Mutation_96Signature.png",
    width=8,height=3,units = "in",res=1200)
sig96_mut
dev.off()

png("~/SNU601/Mutation_StrandBias.png",
    width=8,height=6,units = "in",res=1200)
strandbias
dev.off()

png("~/SNU601/Mutation_GenomeLoc.png",
    width=10,height=4,units = "in",res=1200)
plot_rainfall(mutgr,title = "scATAC-SNU601",chromosomes = chromosomes, 
              cex = 1.5, ylim = 1e+09)
dev.off()

png("~/SNU601/Mutation_GenomeDistribution.png",width=10,height=4,
    units = "in",res=1200)
myplot_spectrum_region(type_occurrences_region)
dev.off()
