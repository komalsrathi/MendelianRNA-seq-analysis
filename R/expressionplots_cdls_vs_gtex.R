# Author: Komal S. Rathi
# Date: 05/13/2019
# Function: Barplots/Boxplots between 4 groups (GTEx whole blood, EBV and CdLS patients positives and negatives)

library(RDiseaseXpress)
library(dplyr)
library(ggplot2)
library(ggpubr)

setwd('~/Projects/DGD_Mendelian_RNASeq/')
source('R/pubTheme.R')
negative <- 'CDL-219-P'

genes <- read.delim('data/CDL_genes.txt', stringsAsFactors = F, header = F)
genes <- genes$V1
if(file.exists('data/expression_data/GTEX_diseaseXpress_expression.RData')){
  print("Exists")
  load('data/expression_data/GTEX_diseaseXpress_expression.RData')
} else {
  dat <- getDataAnnotationByGeneSymbol(myGeneSymbols = genes, myStudy = "GTEx", myNorms = "rsem")
  save(dat, file = 'data/expression_data/GTEX_diseaseXpress_expression.RData')
}
dat <- unique(dat[,c("gene_symbol","gene_id","data.sample_id","data.rsem.fpkm")])
annot <- unique(dat[,c('gene_id','gene_symbol')])

# gtex
setwd('~/Projects/DGD_Mendelian_RNASeq/')
meta <- read.delim('data/gtex_data/GTEx_metadata_593_RNAseq.txt')
dat.sub <- merge(dat, meta, by.x = 'data.sample_id', by.y = 'Run')
dat.sub <- dat.sub[,c('gene_id','data.rsem.fpkm','body_site')]
dat.sub <- dat.sub %>% 
  group_by(gene_id, body_site) %>%
  mutate(mean = mean(data.rsem.fpkm), FPKM = data.rsem.fpkm, Sample = body_site) %>%
  as.data.frame()
dat.sub$Sample <- as.character(dat.sub$Sample)
dat.sub$Sample[dat.sub$Sample == "Cells - EBV-transformed lymphocytes"] <- "EBV-transformed"
dat.sub$color <- "GTEx"
for.bar <- unique(dat.sub[,c('gene_id','mean','Sample','color')])
colnames(for.bar) <- c("gene_id","FPKM",'Sample','color')
for.box <- dat.sub[,c('gene_id','FPKM','Sample','color')]

# patient
p <- readRDS('results/expression/CDLS_patient_expr_hg38.RDS')
p <- p[which(p$gene_id %in% dat.sub$gene_id),c('gene_id','FPKM','Sample')]
p <- p[which(p$Sample != "CDL-219-P"),]
negs <- c('CDL-022-P','CDL-069-99P','CDL-086-99P')
p$color <- ifelse(p$Sample %in% negs,'negative','positive')

# expression barplot
# for each gene, create expression plot across all samples
# color by negative and positive samples
for.bar <- rbind(for.bar, p)
for.bar <- merge(for.bar, annot, by.x = 'gene_id')
for.bar$color[for.bar$color == "positive"] <- "CdLS (positive)"
for.bar$color[for.bar$color == "negative"] <- "CdLS (negative)"
for.bar$color <- ifelse(for.bar$color == "GTEx", for.bar$Sample, for.bar$color)
for.bar$color <- factor(for.bar$color, levels = c("CdLS (positive)","CdLS (negative)","EBV-transformed","Whole Blood"))
bar.plot <- ggplot(data = for.bar, aes(x = Sample, y = FPKM, fill = color)) + 
  geom_bar(stat = 'identity') +
  facet_wrap(~gene_symbol, scales = 'free', ncol = 5) + 
  theme_bw() + theme_Publication() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
bar.plot  
ggsave(bar.plot, filename = 'paper/expression_analysis/CDL_genes_expression_barplot.pdf', device = 'pdf', width = 20, height = 12)  

# boxplot comparing expression in 4 groups
for.box <- rbind(for.box, p)
for.box <- merge(for.box, annot, by.x = 'gene_id')
for.box$Sample_Type <- for.box$color
for.box$Sample_Type <- ifelse(for.box$Sample_Type == "GTEx", for.box$Sample, for.box$Sample_Type)
for.box$Sample_Type[for.box$Sample_Type == "positive"] <- "CdLS (positive)"
for.box$Sample_Type[for.box$Sample_Type == "negative"] <- "CdLS (negative)"
for.box$Sample_Type <- factor(for.box$Sample_Type, levels = c("CdLS (positive)","CdLS (negative)","EBV-transformed","Whole Blood"))
# ggplot(for.box, aes(x = Sample_Type, y = FPKM, fill = Sample_Type)) + 
#   stat_boxplot(geom ='errorbar', width = 0.2) +
#   geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
#   geom_jitter(width = 0.1, pch = 21) +
#   theme_Publication() +
#   theme(axis.text.x = element_blank()) +
#   facet_wrap(~gene_symbol, ncol = 5, scales = "free") +
#   ggtitle("Gene Expression: CDLS Patients (positive and negative) vs \nGTEx Whole Blood and EBV transformed lymphocytes") +
#   stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", label.y = 30) +
#   stat_compare_means(method = "anova", label.y = 20)

my_comparisons <- list(c("CdLS (negative)","CdLS (positive)"),
                       c("CdLS (negative)","EBV-transformed"),
                       c("CdLS (negative)","Whole Blood"))
p <- ggplot(for.box, aes(x = Sample_Type, y = FPKM, fill = Sample_Type)) + 
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
  geom_jitter(width = 0.1, pch = 21) +
  theme_Publication() +
  theme(axis.text.x = element_blank()) +
  facet_wrap(~gene_symbol, ncol = 5, scales = "free") +
  ggtitle("Gene Expression: CDLS Patients (positive and negative) vs \nGTEx Whole Blood and EBV transformed lymphocytes") +
  stat_compare_means(comparisons = my_comparisons)
p
ggsave(p, filename = 'results/expression/CDL_genes_expression_boxplot_4groups.pdf', device = 'pdf', width = 17, height = 15)  
