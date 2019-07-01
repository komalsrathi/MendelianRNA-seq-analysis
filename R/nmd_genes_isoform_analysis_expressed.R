# Author: Komal S. Rathi
# Date: 04/15/2019
# Function: 
# NMD Genes isoform analysis + annotate with omim phenotypic series
# expressed genes only

library(VennDiagram)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)

setwd('~/Projects/DGD_Mendelian_RNASeq/')
load('data/mend_genes_isoform_analysis_input.RData')
source('R/pubTheme.R')

# annot
annot <- read.delim('~/Projects/toil-rnaseq-20k/data/annotation/gencode.v23.annotation_gi_ti_gs.txt')
annot <- unique(annot[,c(1,2)])

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

# add known canonical transcript information
dat <- read.delim('data/knownCanonicalV23.txt')
dat <- dat[which(dat$gene_id %in% patient.mend.genes$gene_id),]
dat$canonical_transcript <- 'Yes'
res <- merge(tpm.data, dat, by = c('gene_id','transcript_id'), all.x = TRUE)
res$canonical_transcript[is.na(res$canonical_transcript)] <- "No"
res <- merge(annot, res, by = 'gene_id')
table <- unique(res[,c('gene_symbol','gene_id','transcript_id','canonical_transcript','expressed_transcript','meanTPM','meanCov')])
table <- table[order(table$gene_id, table$canonical_transcript, decreasing = T),]
write.table(table, file = 'paper/mendelian_disorders_genes/nmd_genes_isoform_expression_CdLS_patients.txt', quote = F, sep = "\t", row.names = F)
table$ensembl <- gsub('[.].*','',table$gene_id)

# barplot of canonical vs alternate splice isoform expression
barchart <- unique(table[,c('gene_id','canonical_transcript','expressed_transcript')])
barchart <- dcast(barchart, gene_id~canonical_transcript, value.var = 'expressed_transcript')
barchart$type <- NA
barchart$type[which(barchart$N == 0 & barchart$Y > 0)] <- 'Canonical'
barchart$type[which(barchart$N > 0 & barchart$Y == 0)] <- 'Non-Canonical'
barchart$type[is.na(barchart$type)] <- 'Both'
barchart <- plyr::count(barchart$type)

p <- ggplot(barchart, aes(x = x, y = freq, fill = x)) + 
  geom_bar(stat = 'identity', color = 'black',  width=0.8) +
  theme_bw() + 
  theme_Publication2() + ylab('# of Expressed Mendelian Genes') + xlab('') +
  ggtitle('Canonical Isoform Expression (NMD Genes)\nFPKM > 1 & Coverage > 10') + 
  geom_text(aes(label = freq, vjust = -0.5, hjust = 0.5), position = position_dodge(width = 0.8), size = 3) +
  guides(fill = FALSE)
p
ggsave(filename = 'paper/mendelian_disorders_genes/nmd_genes_canonical_isoform_expression.png', plot = p, device = 'png', width = 6, height = 4)

# add phenotype data to table
add.pheno <- read.delim('data/genemap_phenotypicSeries_merged.txt', stringsAsFactors = F, check.names = FALSE)
table <- merge(table, add.pheno, by.x = 'ensembl', by.y = 'Ensembl Gene ID', all.x = TRUE)
table <- table %>% group_by(Phenotypic_Superset) %>% 
  mutate(ngenes = length(unique(gene_symbol)),
         perc = ngenes/total.genes*100) %>%
  as.data.frame()
table <- unique(table[,c('gene_symbol','gene_id','MIM Number','Phenotypic Series Number','Phenotypic_Superset','Phenotypes','Mim Number','total.genes','ngenes','perc')])
table <- table[order(table$Phenotypic_Superset, table$perc, decreasing = T),]
colnames(table)[3] <- 'Phenotype MIM number'
colnames(table)[7] <- 'Gene/Locus MIM number'
write.table(table, file = 'paper/mendelian_disorders_genes/nmd_genes_isoforms_with_omim_phenotype_information_1745Genes.txt', quote = F, sep = "\t", row.names = F)
