# Author: Komal S. Rathi
# Date: 04/24/2018
# Function: Create meta file for GTEx RNA-seq samples

setwd('~/Projects/DGD_Mendelian_RNASeq/')
tissues <- c("Cells - EBV-transformed lymphocytes", 
             "Whole Blood")

# sra run table
dat <- read.delim('~/Projects/DGD_Mendelian_RNASeq/data/gtex_data/SraRunTable.txt')
dat <- dat[which(dat$body_site_s %in% tissues),]

# get tissues from GTEx
srr <- read.delim('~/Projects/DGD_Mendelian_RNASeq/data/gtex_data/GTEx_samplesheet_v6.txt')
srr <- srr[which(srr$SMTSD %in% tissues),]

# total RNA-seq files from project
meta <- read.delim('~/Projects/DGD_Mendelian_RNASeq/data/gtex_data/GTEx_1825.txt')
meta <- meta[which(meta$Assay_Type == "RNA-Seq"),]
meta <- meta[,c('Run','Sample_Name','body_site','histological_type')]

setdiff(srr$SAMPID, meta$Sample_Name)
setdiff(dat$Run_s, meta$Run)

# couldn't find these: 
# SRR660957 
# SRR816039
# use 593 samples
write.table(meta, file = '~/Projects/DGD_Mendelian_RNASeq/data/gtex_data/GTEx_metadata_593_RNAseq.txt', quote = F, sep = "\t", row.names = F)


