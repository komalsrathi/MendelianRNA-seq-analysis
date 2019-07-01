# Author: Komal S. Rathi
# Date: 04/15/2019
# Function: 
# Mendelian Genes isoform analysis on all mendelian genes (irrespective of expression)

library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)

# repeat with all 2400 mendelian genes
mend.genes <- read.delim('data/Mendelian_Disorders_genelist_geneids.txt', stringsAsFactors = F, header = F)
colnames(mend.genes)[1] <- 'gene_id'
mend.genes <- mend.genes[-which(mend.genes$gene_id %in% c('ENSG00000274276.4','ENSG00000255330.8')),]

# annotation
annot <- read.delim('~/Projects/toil-rnaseq-20k/data/annotation/gencode.v23.annotation_gi_ti_gs.txt')
annot <- unique(annot[,c(1,2)])
annot <- annot[which(annot$gene_id %in% mend.genes),]

# cdls data
load('data/mend_genes_isoform_analysis_input.RData')
patient <- patient[which(patient$mend.gene == "Y"),c('gene_symbol','type')]
gtex.ebv <- gtex.ebv[which(gtex.ebv$mend.gene == "Y"),c('gene_symbol','type')]
annot$EBV_expr <- ifelse(annot$gene_symbol %in% gtex.ebv$gene_symbol,'Y','N')
annot$CdLS_expr <- ifelse(annot$gene_symbol %in% patient$gene_symbol,'Y','N')

# take the max of the mean TPM to resolve duplicates
tpm.data <- readRDS('data/expression_data/CDLS_patient_expr_transcripts_hg38.RDS')
tpm.data <- tpm.data[which(tpm.data$Sample != 'CDL-219-P' & tpm.data$gene_id %in% mend.genes),]
tpm.data <- merge(annot, tpm.data, by = 'gene_id')
tpm.data <- tpm.data %>% group_by(transcript_id, gene_id, gene_symbol) %>%
  summarise(meanTPM = mean(TPM),
            meanCov = mean(expected_count),
            expressed_transcript = ifelse(meanCov > 10 & meanTPM > 1, "Yes", "No")) %>%
  as.data.frame()
test <- tpm.data %>% group_by(gene_symbol, gene_id) %>% summarise(meanTPM = mean(meanTPM))
test <- test[order(test$gene_symbol, test$meanTPM, decreasing = T),]
test <- test[!duplicated(test$gene_symbol),]
annot <- annot[which(annot$gene_id %in% test$gene_id),]

# now add the annotation to this table
annot$ensembl <- gsub('[.].*','',annot$gene_id)
add.pheno <- read.delim('data/genemap_phenotypicSeries_merged.txt', stringsAsFactors = F)
table <- merge(annot, add.pheno, by.x = 'ensembl', by.y = 'Ensembl.Gene.ID', all.x = TRUE)
EBV.expr <- table %>% group_by(Phenotypic_Superset) %>%
  filter(EBV_expr == "Y") %>%
  mutate(EBV.genes = length(unique(gene_symbol)),
         EBV.perc = EBV.genes/total.genes*100) %>%
  select(Phenotypic_Superset, EBV.genes, EBV.perc) %>%
  unique() %>%
  as.data.frame()
CdLS.expr <- table %>% group_by(Phenotypic_Superset) %>%
  filter(CdLS_expr == "Y") %>%
  mutate(CdLS.genes = length(unique(gene_symbol)),
         CdLS.perc = CdLS.genes/total.genes*100) %>%
  select(Phenotypic_Superset, CdLS.genes, CdLS.perc) %>%
  unique() %>%
  as.data.frame()
tmp <- merge(EBV.expr, CdLS.expr, by = 'Phenotypic_Superset', all = TRUE)
final <- merge(table, tmp, by = 'Phenotypic_Superset', all.x = TRUE)
colnames(final)[c(7,10)] <- c('Phenotype MIM number','Gene/Locus MIM number')
final$ensembl <- NULL
final <- final[,c('Phenotypic_Superset','Phenotype MIM number','Phenotypic.Series.Number','Phenotypes','Gene/Locus MIM number','gene_id','gene_symbol','total.genes','EBV_expr','EBV.genes','EBV.perc','CdLS_expr','CdLS.genes','CdLS.perc')]
write.table(final, file = 'paper/mendelian_disorders_genes/nmd_genes_isoforms_with_omim_phenotype_information_2541Genes.txt', quote = F, sep = "\t", row.names = F)
