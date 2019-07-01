# Author: Komal S. Rathi
# Date: 04/15/2019
# Function: 
# For all 1745 expressed genes in CdLS
# Isoform comparisons between Brain, EBV and CdLS

library(VennDiagram)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)

setwd('~/Projects/DGD_Mendelian_RNASeq/')
load('data/mend_genes_isoform_analysis_input.RData')
source('R/pubTheme.R')

# isoform expression table and barplot
# for 1745 genes expressed in patient, get TPM values
patient.mend.genes <- merge(patient, p.expr, by = 'gene_symbol')
patient.mend.genes <- patient.mend.genes[which(patient.mend.genes$mend.gene == "Y"),]

# isoform data
# 6460/18385 transcripts are expressed
tpm.data <- readRDS('data/expression_data/CDLS_patient_expr_transcripts_hg38.RDS')
tpm.data <- tpm.data[which(tpm.data$Sample != 'CDL-219-P' & tpm.data$gene_id %in% patient.mend.genes$gene_id),]
tpm.data <- tpm.data %>% group_by(transcript_id, gene_id) %>%
  mutate(meanTPM = mean(TPM),
         meanCov = mean(expected_count),
         expressed_transcript = ifelse(meanCov > 10 & meanTPM > 1, "Yes", "No")) %>%
  as.data.frame()

# test 1 at least 1 transcript in each gene in expressed
test <- unique(tpm.data[,c('gene_id','expressed_transcript')]) 
if(length(unique(test$gene_id)) == 1745){
  print("Test passed")
}

# correlation of max transcript expression per gene
# Brain, EBV and CdLS
cdls <- unique(tpm.data[,c('gene_id','transcript_id','meanCov','meanTPM')])
colnames(cdls)[c(3,4)] <- paste0('CdLS_',colnames(cdls)[c(3,4)])
brain <- readRDS('data/Brain_meanTPMandCoverage.RDS')
brain <- brain[which(brain$gene_id %in% patient.mend.genes$gene_id),]
colnames(brain)[c(3,4)] <- paste0('Brain_',colnames(brain)[c(3,4)])
ebv <- readRDS('data/EBV_meanTPMandCoverage.RDS')
ebv <- ebv[which(ebv$gene_id %in% patient.mend.genes$gene_id),]
ebv$body_site <- NULL
colnames(ebv)[c(3,4)] <- paste0('EBV_',colnames(ebv)[c(3,4)])

tpm.comp <- merge(ebv, brain, by = c('gene_id','transcript_id','gene_symbol'))
tpm.comp <- merge(tpm.comp, cdls, by = c('gene_id','transcript_id'))
tpm.comp <- tpm.comp %>% group_by(gene_id) %>%
  mutate(EBV_max = ifelse(EBV_meanTPM == max(EBV_meanTPM), 1, 0),
         Brain_max = ifelse(Brain_meanTPM == max(Brain_meanTPM), 1, 0),
         CdLS_max = ifelse(CdLS_meanTPM == max(CdLS_meanTPM), 1, 0)) %>%
  as.data.frame()
nrow(tpm.comp[which(tpm.comp$EBV_max == 1 & tpm.comp$Brain_max == 1),]) # 1134/1745
nrow(tpm.comp[which(tpm.comp$CdLS_max == 1 & tpm.comp$Brain_max == 1),]) # 1082/1745
nrow(tpm.comp[which(tpm.comp$EBV_max == 1 & tpm.comp$CdLS_max == 1),]) # 1568/1745

# venn diagram of max expressed transcript
png(filename = 'paper/mendelian_disorders_genes/nmd_genes_max_expressed_transcript.png', width = 6, height = 5, units = "in", res = 240)
one <- tpm.comp[which(tpm.comp$CdLS_max == 1),'transcript_id']
two <- tpm.comp[which(tpm.comp$EBV_max == 1),'transcript_id']
three <- tpm.comp[which(tpm.comp$Brain_max == 1),'transcript_id']
p <- draw.triple.venn(area1 = length(one), 
                      area2 = length(two), 
                      area3 = length(three), 
                      n12 = length(intersect(one, two)),
                      n23 = length(intersect(two, three)), 
                      n13 = length(intersect(one, three)),
                      n123 = length(intersect(intersect(one, two), three)), 
                      scaled = F, euler.d = F,
                      category = c(paste0("CdLS Patients\n(",length(one),")"), 
                                   paste0("EBV\n(",length(two),")"), 
                                   paste0("Brain\n(",length(three),")")),
                      fill = c("skyblue", "pink1", "mediumorchid"), lty = "blank")
grid.arrange(gTree(children=p), top="NMD Genes\nMax Expressed Transcript")
dev.off()
# ggsave(filename = 'paper/mendelian_disorders_genes/max_expressed_transcript.png', device = 'png', plot = p, width = 7, height = 5)

# venn diagram of all ~18k transcripts
tpm.comp <- tpm.comp %>% group_by(gene_id, transcript_id) %>%
  mutate(EBV_expr = ifelse(EBV_meanTPM > 1 & EBV_meanCov > 10, 1, 0),
         Brain_expr = ifelse(Brain_meanTPM > 1 & Brain_meanCov > 10, 1, 0),
         CdLS_expr = ifelse(CdLS_meanTPM > 1 & CdLS_meanCov > 10, 1, 0)) %>%
  as.data.frame()
write.table(tpm.comp, file = 'paper/mendelian_disorders_genes/nmd_genes_isoform_expression_1745genes.txt', quote = F, sep = "\t", row.names = F)

png(filename = 'paper/mendelian_disorders_genes/nmd_genes_all_expressed_transcripts.png', width = 6, height = 5, units = "in", res = 240)
one <- tpm.comp[which(tpm.comp$CdLS_expr == 1),'transcript_id']
two <- tpm.comp[which(tpm.comp$EBV_expr == 1),'transcript_id']
three <- tpm.comp[which(tpm.comp$Brain_expr == 1),'transcript_id']
p <- draw.triple.venn(area1 = length(one), 
                      area2 = length(two), 
                      area3 = length(three), 
                      n12 = length(intersect(one, two)),
                      n23 = length(intersect(two, three)), 
                      n13 = length(intersect(one, three)),
                      n123 = length(intersect(intersect(one, two), three)), 
                      scaled = F, euler.d = F,
                      category = c(paste0("CdLS Patients\n(",length(one),")"), 
                                   paste0("EBV\n(",length(two),")"), 
                                   paste0("Brain\n(",length(three),")")),
                      fill = c("skyblue", "pink1", "mediumorchid"), lty = "blank")
grid.arrange(gTree(children=p), top="NMD Genes\nAll Expressed Transcripts")
dev.off()
# ggsave(filename = 'paper/mendelian_disorders_genes/expressed_transcripts.png', device = 'png', plot = p, width = 7, height = 5)

# scatter plot of 18385 transcripts (add vertical and horizontal lines)
# sc.plot <- melt(cov.comp)
# sc.plot <- sc.plot[which(sc.plot$variable %in% c('EBV_TPM','Brain_TPM','CdLS_TPM')),]
p.cor <- cor.test(log2(tpm.comp$EBV_meanTPM+1), log2(tpm.comp$Brain_meanTPM+1), method = "pearson")
p <- ggplot(tpm.comp, aes(log2(EBV_meanTPM+1), log2(Brain_meanTPM+1))) + 
  geom_point(pch = 21, size = 1) + geom_smooth(method=lm, size = 0.5, color = 'darkgray') +
  geom_hline(yintercept = log2(1+1), linetype = "dashed", color = "red") +
  geom_vline(xintercept = log2(1+1), linetype = "dashed", color = "red") + theme_Publication4() +
  xlab('EBV log2(TPM)') + ylab('Brain log2(TPM)') + 
  annotate(geom = "text", x = 5, y = 10, label = paste0("Pearson: 0.74\n Pval:<2.2e-16"), color = "red", size = 3)
p
q.cor <- cor.test(log2(tpm.comp$CdLS_meanTPM+1), log2(tpm.comp$Brain_meanTPM+1), method = "pearson")
q <- ggplot(tpm.comp, aes(log2(CdLS_meanTPM+1), log2(Brain_meanTPM+1))) + 
  geom_point(pch = 21, size = 1) + geom_smooth(method=lm, size = 0.5, color = 'darkgray') +
  geom_hline(yintercept = log2(1+1), linetype = "dashed", color = "red") +
  geom_vline(xintercept = log2(1+1), linetype = "dashed", color = "red") + theme_Publication4() +
  xlab('CdLS log2(TPM)') + ylab('Brain log2(TPM)') +
  annotate(geom = "text", x = 5, y = 10, label = paste0("Pearson: 0.68\n Pval:<2.2e-16"), color = "red", size = 3)
q
r.cor <- cor.test(log2(tpm.comp$CdLS_meanTPM+1), log2(tpm.comp$EBV_meanTPM+1), method = 'pearson')
r <- ggplot(tpm.comp, aes(log2(CdLS_meanTPM+1), log2(EBV_meanTPM+1))) + 
  geom_point(pch = 21, size = 1) + geom_smooth(method=lm, size = 0.5, color = 'darkgray') +
  geom_hline(yintercept = log2(1+1), linetype = "dashed", color = "red") +
  geom_vline(xintercept = log2(1+1), linetype = "dashed", color = "red") + theme_Publication4() +
  xlab('CdLS log2(TPM)') + ylab('EBV log2(TPM)') +
  annotate(geom = "text", x = 5, y = 10, label = paste0("Pearson: 0.96\n Pval:<2.2e-16"), color = "red", size = 3)
r

g <- gridExtra::arrangeGrob(p, q, r, nrow = 1, 
                            top = textGrob('NMD Genes: Correlation of Transcript Expression\n(nGenes = 1745, nTranscripts = 18385)',
                                           gp = gpar(fontsize = 12)))
ggsave(filename = 'paper/mendelian_disorders_genes/nmd_genes_isoform_expression_scatterplots.pdf', device = "pdf", plot = g, width = 10, height = 4)
