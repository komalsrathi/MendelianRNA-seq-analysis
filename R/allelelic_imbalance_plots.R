# Author: Komal S. Rathi
# Date: 06/19/2019
# Function: 
# Boxplots for each gene across all patients to show Allele Imbalance using gnomAD allele balance
# This uses the RData object generated from allelic_imbalance_gatk3.R

library(GenomicRanges)
library(ggplot2)
setwd('~/Projects/DGD_Mendelian_RNASeq/')
source('R/pubTheme.R')

# add AB from gnomad 
if(file.exists('results/allele_imbalance/gnomad_median_AB.RData')){
  load('results/allele_imbalance/gnomad_median_AB.RData')
} else {
  dat <- data.table::fread('results/allele_imbalance/gnomad.txt') # only PASS-ed variants
  dat <- dat[which(dat$AF_POPMAX < 0.002),] # filter by AF_POPMAX
  dat <- as.data.frame(dat)
  dat$id <- paste0('id_',rownames(dat))
  
  # gene list
  cdls <- read.delim('data/CDL_genes.list')
  cdls$score <- 0
  cdls <- unique(cdls[,c('chromosome_name','start_position','end_position','hgnc_symbol','score','strand')])
  
  # overlap using genomic ranges
  subject <- with(cdls, GRanges(chromosome_name, IRanges(start = start_position, end = end_position, names = hgnc_symbol)))
  query <- with(dat, GRanges(CHROM, IRanges(start = POS, end = POS, names = id)))
  
  # find overlaps and subset maf 
  res <- findOverlaps(query = query, subject = subject, type = "within")
  res.df <- data.frame(dat[queryHits(res),], cdls[subjectHits(res),])
  res.df <- unique(res.df[,c('hgnc_symbol','AB_MEDIAN')])
  res.df$label <- 'gnomAD_AB_median'
  gnomad <- res.df
  save(gnomad, file = 'results/allele_imbalance/gnomad_median_AB.RData')
}
colnames(gnomad) <- c("Hugo_Symbol","AB_ref","label")
gnomad$AB_ref <- as.numeric(gnomad$AB_ref)

# allele imbalance from patient data
load('results/allele_imbalance/gatk3_allelic_imbalance.RData')
gatk3$label <- paste0(gatk3$Tumor_Sample_Barcode, " (n = ",gatk3$var.count,")")
gatk3 <- gatk3[,c("Hugo_Symbol","AB_ref","label")]

# total
total <- rbind(gatk3, gnomad)

# genes
genes <- unique(gatk3$Hugo_Symbol)

for(i in 1:length(genes)){
  dat <- total[which(total$Hugo_Symbol == genes[i]),]
  fname <- paste0('results/allele_imbalance/',genes[i],'.png')
  p <- ggplot(dat, aes(x = label, y = AB_ref)) + 
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
    geom_jitter(width = 0.1, pch = 21) +
    # geom_boxplot() + 
    facet_grid(~Hugo_Symbol) + 
    theme_bw() + theme_Publication2() +  
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + ylab('Allele Balance') +
    geom_hline(yintercept = 0.5, color = "red", linetype = "dashed")
  ggsave(p, filename = fname, device = "png", width = 8, height = 6)
}
