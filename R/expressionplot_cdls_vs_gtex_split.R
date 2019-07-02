# Author: Komal S. Rathi
# Date: 04/04/2019
# Function: 
# Boxplots between 3 groups (GTEx whole blood, EBV and CdLS patients)
# Boxplot between CdLS patients with NIPBL splice mutation, without splice mutation and EBV

library(RDiseaseXpress)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd('~/Projects/DGD_Mendelian_RNASeq/')
source('R/pubTheme.R')
negative <- 'CDL-219-P'

# gene expression data
genes <- read.delim('data/CDL_genes.txt', stringsAsFactors = F, header = F)
genes <- genes$V1
if(file.exists('data/expression_data/GTEX_diseaseXpress_expression.RData')){
  load('data/expression_data/GTEX_diseaseXpress_expression.RData')
} else {
  dat <- getDataAnnotationByGeneSymbol(myGeneSymbols = genes, myStudy = "GTEx", myNorms = "rsem")
  save(dat, file = 'data/expression_data/GTEX_diseaseXpress_expression.RData')
}
dat <- unique(dat[,c("gene_symbol","gene_id","data.sample_id","data.rsem.fpkm")])
annot <- unique(dat[,c('gene_id','gene_symbol')])

# meta
setwd('~/Projects/DGD_Mendelian_RNASeq/')
meta <- read.delim('data/gtex_data/GTEx_metadata_593_RNAseq.txt')
dat.sub <- merge(dat, meta, by.x = 'data.sample_id', by.y = 'Run')
dat.sub <- dat.sub[,c('gene_id','data.rsem.fpkm','body_site','data.sample_id')]
colnames(dat.sub) <- c("gene_id","FPKM",'Sample','id')

# barplot
p <- readRDS('results/expression/CDLS_patient_expr_hg38.RDS')
p <- p[which(p$Sample != negative),]
p <- p[which(p$gene_id %in% dat.sub$gene_id),c('gene_id','FPKM','Sample')]
negs <- c('95-0614-P',
          'CDL-022-P',
          'CDL-086-99P',
          'CDL-069-99P',
          'CDL-679-14P')
colnames(p)[3] <- 'id'
p$Sample <- ifelse(p$id %in% negs,'CdLS_negative','CdLS_positive')
total <- rbind(dat.sub, p)
total$Sample <- as.character(total$Sample)
total$Sample[total$Sample == "Cells - EBV-transformed lymphocytes"] <- "EBV-transformed"
total <- merge(total, annot, by.x = 'gene_id')
colnames(total)[3] <- "Groups"

# does NIPBL splice mutation affect its expression? No.
nibpl.splice <- total[which(total$gene_symbol == "NIPBL"),]
nibpl.splice$Groups[nibpl.splice$id %in% c('27571P', '95-1682-P', 'CDL-075-99P', 'CDL-123-01P', 'CDL-217-05P', 'CDL-515-09P')] <- "NIPBL Splice (+)"
nibpl.splice$Groups[nibpl.splice$Groups %in% c("CdLS_positive", "CdLS_negative")] <- 'NIPBL Splice (-)'
nibpl.splice$Groups <- factor(nibpl.splice$Groups, levels = c("NIPBL Splice (+)","NIPBL Splice (-)","EBV-transformed","Whole Blood"))
my_comparisons <- list(c("NIPBL Splice (+)", "EBV-transformed"),
                       c("NIPBL Splice (+)", "Whole Blood"),
                       c("NIPBL Splice (+)", "NIPBL Splice (-)"))
p <- ggplot(nibpl.splice, aes(x = Groups, y = FPKM, fill = Groups)) + 
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, pch = 21) +
  theme_Publication() +
  theme(axis.text.x = element_blank()) +
  ggtitle("NIPBL Gene Expression") +
  stat_compare_means(comparisons = my_comparisons) 
p

# expression boxplot
# split into positive and negative
# colnames(total)[3] <- 'Groups'
pos <- total[which(total$Groups != "CdLS_negative"),]
neg <- total[which(total$Groups != "CdLS_positive"),]

# extra plots
my_comparisons <- list(c("CdLS_positive", "EBV-transformed"),
                       c("CdLS_positive", "Whole Blood"))
q <- ggplot(pos, aes(x = Groups, y = FPKM, fill = Groups)) + 
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, pch = 21) +
  theme_Publication() +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~gene_symbol, ncol = 5, scales = "free") +
  ggtitle("Gene Expression: CDLS Patients vs GTEx Whole Blood and EBV transformed lymphocytes") +
  stat_compare_means(comparisons = my_comparisons) 
q
# ggsave(q, filename = 'paper/expression_analysis/CdLS_positives_vs_GTEx_v2.pdf', device = 'pdf', width = 17, height = 14)  

q <- ggplot(pos, aes(x = Groups, y = FPKM, fill = Groups)) + 
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
  theme_Publication() +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~gene_symbol, ncol = 5, scales = "free") +
  ggtitle("Gene Expression: CDLS Patients vs GTEx Whole Blood and EBV transformed lymphocytes") +
  stat_compare_means(comparisons = my_comparisons) 
q
# ggsave(q, filename = 'paper/expression_analysis/CdLS_positives_vs_GTEx_v3.pdf', device = 'pdf', width = 17, height = 14)  

# separated boxplots
# dodge boxplot (positives)
pd = position_dodge(width = 0.7)
q <- ggplot(pos, aes(x = gene_symbol, y = FPKM, fill = Groups)) + 
  stat_boxplot(geom ='errorbar', position = pd, width = 0.2) +
  geom_boxplot(position = position_dodge(width=0.7), lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
  theme_Publication2() +
  xlab("") +
  ggtitle("Gene Expression: CDLS Patients (+) vs GTEx Whole Blood and EBV transformed lymphocytes") 
q
ggsave(q, filename = 'paper/expression_analysis/CdLS_positives_vs_GTEx.pdf', device = 'pdf', width = 14, height = 6)  

# dodge boxplot (negatives)
pd = position_dodge(width = 0.7)
q <- ggplot(neg, aes(x = gene_symbol, y = FPKM, fill = Groups)) + 
  stat_boxplot(geom ='errorbar', position = pd, width = 0.2) +
  geom_boxplot(position = position_dodge(width=0.7), lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
  theme_Publication2() +
  xlab("") +
  ggtitle("Gene Expression: CDLS Patients (-) vs GTEx Whole Blood and EBV transformed lymphocytes") 
q
ggsave(q, filename = 'paper/expression_analysis/CdLS_negatives_vs_GTEx.pdf', device = 'pdf', width = 14, height = 6)  


# transcript level changes
# p <- readRDS('results/expression/CDLS_patient_expr_transcripts_hg38.RDS')
# p <- p[which(p$gene_id %in% annot$gene_id),]
# 
# for(i in 1:nrow(annot)){
#   gi <- annot[i,1]
#   gs <- annot[i,2]
#   tmp <- p[which(p$gene_id %in% gi),]
#   fname <- paste0('results/expression/',gs,'_transcripts.pdf')
#   t <- ggplot(tmp, aes(transcript_id, TPM)) + 
#     geom_bar(stat = 'identity') + theme_Publication() +
#     theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
#     facet_wrap(~Sample, ncol = 4, scales = "free") +
#     ggtitle(paste0(gs," (TPM)")) + ylab("") +
#     theme(plot.margin = unit(c(.5,.5,.5,2), "cm"))
#   if(i == 3){
#     ggsave(filename = fname, plot = t, device = 'pdf', width = 25, height = 10)
#   } else {
#     ggsave(filename = fname, plot = t, device = 'pdf', width = 15, height = 10)
#   }
# }

# ggplot(pos, aes(x = Groups, y = FPKM, fill = Groups)) + 
#   stat_boxplot(geom ='errorbar', width = 0.2) +
#   geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
#   theme_Publication5() +
#   theme(axis.text.x = element_blank()) + xlab('') +
#   facet_wrap(strip.position = "bottom",~gene_symbol, nrow = 1) +
#   ggtitle("Gene Expression: CDLS Patients vs GTEx Whole Blood and EBV transformed lymphocytes") +
#   stat_compare_means(comparisons = my_comparisons, size = 3, tip.length = 0.01)
