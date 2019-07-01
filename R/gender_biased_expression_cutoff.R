# Author: Komal S. Rathi
# Date: 04/04/2019
# Function: Determine the FPKM cut-off value to be used to determine expression
# using only protein coding gender specific genes

library(dplyr)
library(reshape2)
library(ggplot2)
setwd('~/Projects/DGD_Mendelian_RNASeq/')
source('R/pubTheme.R')
load('data/expression_data/gtex_cdls_matrix.RData')

# annotation and get all Y-specific genes
annot <- data.table::fread('~/Projects/toil-rnaseq-20k/data/annotation/gencode.v23.annotation.txt')
annot <- unique(annot[which(annot$chr %in% c("chrY", "chrX") & annot$biotype == "protein_coding"),c("gene_symbol","chr","biotype")])
annot <- annot %>%
  group_by(gene_symbol) %>%
  dplyr::summarise(chr = toString(chr),
            chr.ct = n()) %>%
  filter(chr.ct == 1 & chr %in% c("chrY","chrX"))

# get patient data
patients <- dat[,grep('^SRR', colnames(dat), invert = T)]

# format for sex check
if(file.exists('data/patient_genes.fpkm.gct')){
  print("File exists")
} else {
  annot <- data.table::fread('~/Projects/toil-rnaseq-20k/data/annotation/gencode.v23.annotation.txt')
  annot <- unique(annot[,c('gene_id','gene_symbol')])
  annot <- as.data.frame(annot)
  for.gender <- merge(annot, patients, by.x = 'gene_symbol', by.y = 'row.names')
  for.gender <- for.gender[!duplicated(for.gender$gene_symbol),]
  for.gender <- for.gender[,c(2,1,3:ncol(for.gender))]
  colnames(for.gender)[1:2] <- c("Name","Description")
  df <- t(data.frame(c("#1.2", rep(NA,16)),
                     c("58581", 15, rep(NA,15)), 
                     colnames(for.gender),
                     stringsAsFactors = F))
  df <- data.frame(df, stringsAsFactors = F)
  rownames(df) <- NULL
  colnames(df) <- colnames(for.gender)
  df <- rbind(df, for.gender)
  df[is.na(df)] <- ''
  write.table(df, file = 'data/patient_genes.fpkm.gct', quote = F, sep = "\t", row.names = F, col.names = F)
}

sex_biased_genes <- read.delim("data/sex_biased_genes.txt",header=T,stringsAsFactors = F,strip.white=T)
f.genes <- sex_biased_genes[which(sex_biased_genes$escape.X.inactivation == "YES" & sex_biased_genes$chr == "chrX"),'gene_name']
m.genes <- sex_biased_genes[which(sex_biased_genes$chr == "chrY"),'gene_name']
f.genes <- f.genes[f.genes %in% annot$gene_symbol] # 16
m.genes <- m.genes[m.genes %in% annot$gene_symbol] # 13

# boxplot of expression of these genes in male patients and female patients
patients <- melt(as.matrix(patients))
genderinfo <- read.delim('data/gender_info.txt')
patients <- merge(patients, genderinfo, by.x = 'Var2', by.y = 'Patient')
patients$Label[patients$Var1 %in% f.genes] <- 'Female_Genes'
patients$Label[patients$Var1 %in% m.genes] <- 'Male_Genes'
patients <- patients[!is.na(patients$Label),]

m <- patients[which(patients$Label == "")]

p <- ggplot(data = patients, 
       aes(x = Gender, y = log2(value+1), color = Gender)) + 
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, pch = 21) +
  facet_wrap(~Label, scales = "free_y") +
  theme_bw() + theme_Publication() +
  ylab('log2(FPKM + 1)') +
  ggtitle('Expression Cut-off') + guides(color = F)
ggsave(p, filename = 'paper/expression_analysis/expr_cutoff.pdf', device = 'pdf', width = 8, height = 6)  

# This shows that you can use FPKM > 1 as expressed cutoff and
# any gene that has 0 FPKM value across all samples can be used as not-expressed