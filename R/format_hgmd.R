# Author: Komal S. Rathi
# Date: 01/31/2019
# Function: 
# convert hgmd vcf to tab delimited file (using hgmd_convert_vcf2tab.sh)
# 1. format the above file to only keep disease causing variants from HGMD (DM or DM?)
# 2. subset to CdLS genes only

setwd('/mnt/isilon/cbmi/variome/rathik/mendelian_rnaseq')
library(data.table)
library(reshape2)
dat <- data.table::fread('hgmd/out.INFO') # delimited file from HGMD vcf
dat <- as.data.frame(dat)
dat <- cbind(dat, colsplit(dat$DNA, pattern = ":", names = c('RefSeq','HGVSc')))
dat <- cbind(dat, colsplit(dat$PROT, pattern = ":", names = c('RefSeq_Protein','HGVSp_Short')))
dat$DNA <- NULL
dat$PROT <- NULL
dat <- dat[grep('DM', dat$CLASS),] # only disease causing variants 

# subset to cdls genes
cdls <- read.delim('MendelianRNA-seq/data/CDL_genes.txt', stringsAsFactors = F, header = F)
dat <- dat[which(dat$GENE %in% cdls$V1),]
write.table(dat, file = 'hgmd/hgmd_pro_2019.1_hg19.tab', quote = F, sep = "\t", row.names = F)
