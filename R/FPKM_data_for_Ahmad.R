setwd('~/Projects/DGD_Mendelian_RNASeq/')

library(dplyr)
library(reshape2)

load('data/mend_genes_isoform_analysis_input.RData')
rm(list = grep('gtex.ebv|gtex.wb', ls(), invert = T, value = T))
gtex.nmd <- gtex.wb %>%
  filter(mend.gene == "Y")

# annotation
annot <- data.table::fread('~/Projects/toil-rnaseq-20k/data/annotation/gencode.v23.annotation.txt')
annot <- unique(annot[which(annot$biotype == "protein_coding"),c('gene_id','gene_symbol')])
annot <- as.data.frame(annot)
annot$ensembl <- gsub('[.].*','',annot$gene_id)

expr.genes <- function(rds.path){
  dat <- readRDS(rds.path)
  dat <- acast(dat, gene_id~id, value.var = 'expected_count')
  dat <- dat[apply(dat, MARGIN = 1, FUN = function(x) all(x >= 10)),]
  dat <- as.data.frame(dat)
  dat <- merge(annot[,1:2], dat, by.x = 'gene_id', by.y = 'row.names')
  dat$means <- rowMeans(dat[,3:ncol(dat)])
  dat <- dat[order(dat$gene_symbol, dat$means, decreasing = TRUE),]
  dat <- dat[!duplicated(dat$gene_symbol),]
  # expr.genes <- unique(dat$gene_symbol)
  expr.genes <- unique(dat[,c('gene_id','gene_symbol')])
  return(expr.genes)
}

# gene id and symbol
p.expr <- expr.genes(rds.path = 'data/expression_data/CDLS_patient_expr_hg38.RDS')
wb.expr <- expr.genes(rds.path = 'data/expression_data/WB_expr_hg38.RDS')
ebv.expr <- expr.genes(rds.path = 'data/expression_data/EBV_expr_hg38.RDS')

# gtex ebv total
gtex.ebv <- gtex.ebv %>%
  inner_join(annot, by = "gene_symbol") %>%
  filter(gene_id %in% ebv.expr$gene_id)

# gtex wb total
gtex.wb <- gtex.wb %>%
  inner_join(annot, by = "gene_symbol") %>%
  filter(gene_id %in% wb.expr$gene_id)

# gtex wb NMD
gtex.nmd <- gtex.nmd %>%
  inner_join(annot, by = "gene_symbol") %>%
  filter(gene_id %in% wb.expr$gene_id)

# expression
wb.allexpr <- readRDS('data/expression_data/WB_expr_hg38.RDS')
wb.total <- wb.allexpr %>%
  filter(gene_id %in% gtex.wb$gene_id) %>% 
  inner_join(annot[,c('gene_id','gene_symbol')], by = "gene_id") %>%
  select(gene_id, gene_symbol, FPKM, id) %>%
  spread(id, FPKM)
write.table(wb.total, file  = 'results/expression/WB_total_FPKM_matrix.txt', quote = F, sep = "\t", row.names = F)

wb.nmd <- wb.allexpr %>%
  filter(gene_id %in% gtex.nmd$gene_id) %>% 
  inner_join(annot[,c('gene_id','gene_symbol')], by = "gene_id") %>%
  select(gene_id, gene_symbol, FPKM, id) %>%
  spread(id, FPKM)
write.table(wb.nmd, file  = 'results/expression/WB_NMD_FPKM_matrix.txt', quote = F, sep = "\t", row.names = F)

ebv.allexpr <- readRDS('data/expression_data/EBV_expr_hg38.RDS')
ebv.total  <- ebv.allexpr %>%
  filter(gene_id %in% gtex.ebv$gene_id) %>% 
  inner_join(annot[,c('gene_id','gene_symbol')], by = "gene_id") %>%
  select(gene_id, gene_symbol, FPKM, id) %>%
  spread(id, FPKM)
write.table(ebv.total, file  = 'results/expression/EBV_total_FPKM_matrix.txt', quote = F, sep = "\t", row.names = F)

