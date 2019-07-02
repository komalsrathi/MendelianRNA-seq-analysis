# Author: Komal S Rathi
# Date: 01/31/2019
# Function: 
# 1. subset to only CDLS genes
# 2. add HGMD annotation to MAF
# This will generate *-hgmdannotated.maf files 
# which will be used in various scripts

library(data.table)
library(hutils)
library(tidyr)
library(GenomicRanges)
library(reshape2)
library(dplyr)

setwd('/mnt/isilon/cbmi/variome/rathik/mendelian_rnaseq')

# hgmd annotation in tab-delimited format
hgmd <- read.delim('hgmd/hgmd_pro_2019.1_hg19.tab', stringsAsFactors = F)
hgmd <- unique(hgmd[,c('CHROM','POS','REF','ALT','GENE','STRAND','HGVSc','CLASS')])

# CdLS gene list
cdls <- read.delim('MendelianRNA-seq/data/CDL_genes.txt', header = F, stringsAsFactors = F)
cdls <- cdls$V1

# source directory
folder <- 'updated-mafs'

lf <- list.files(path = folder, pattern = '*.maf', full.names = TRUE, recursive = TRUE)
for(i in 1:length(lf)){
  print(lf[i])

  newfilename <- gsub('[.]maf','-hgmdannotated.maf',lf[i])
  print(newfilename)

  dat <- data.table::fread(lf[i], verbose = FALSE)
  dat <- as.data.frame(dat)

  # subset to CdLS genes
  dat <- dat[which(dat$Hugo_Symbol %in% cdls),] 

  print(dim(dat))

  # merge with HGMD
  dat <- merge(x = dat, y = hgmd, 
    by.x = c("Chromosome","Start_Position","Reference_Allele","Hugo_Symbol","Strand","HGVSc"), 
    by.y = c("CHROM","POS","REF","GENE","STRAND","HGVSc"), all.x = TRUE)

  print(dim(dat))

  # write out
  write.table(dat, file = newfilename, quote = F, sep = "\t", row.names = F)
}



