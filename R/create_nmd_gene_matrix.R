# Author: Komal S. Rathi
# Date: 04/04/2019
# Function: Create Heatmap of NMD Genes

library(reshape2)
library(dplyr)
library(pheatmap)
library(ggplot2)

setwd('~/Projects/DGD_Mendelian_RNASeq/')

# matrix of mendelian genes
if(file.exists('data/expression_data/mendelian_genes_matrix.RData')){
  print("File exists")
  load('data/expression_data/mendelian_genes_matrix.RData')
} else {
  p <- readRDS('data/expression_data/CDLS_patient_expr_hg38.RDS')
  p <- p[which(p$id != 'CDL-219-P'),]
  n <- readRDS('data/expression_data/GTEx_Blood_442.RDS')
  
  n.meta <- read.delim('data/expression_data/GTEx_metadata.txt')
  p.meta <- data.frame(id = unique(p$id), type = "CdLS_Patient")
  p.meta <- p.meta[which(p.meta$id != 'CDL-219-P'),]
  meta <- rbind(n.meta, p.meta)
  
  dat <- rbind(p, n)
  counts <- dcast(dat, gene_id~id, value.var = 'expected_count')
  dat <- dcast(dat, gene_id~id, value.var = 'FPKM')
  annot <- data.table::fread('~/Projects/toil-rnaseq-20k/data/annotation/gencode.v23.annotation_gi_ti_gs.txt')
  annot <- annot[,1:2]
  tmp <- plyr::count(annot$gene_symbol)
  tmp <- tmp[which(tmp$freq > 1),]
  
  dat <- merge(annot, dat, by = 'gene_id')
  counts <- merge(annot, counts, by = 'gene_id')
  
  ncol <- ncol(dat)
  dat <- dat %>% 
    as.data.frame() %>%
    #dplyr::select(-gene_id) %>% 
    unique() %>% 
    mutate(means = rowMeans(.[3:ncol])) %>% 
    arrange(desc(means)) %>% 
    distinct(gene_symbol, .keep_all = TRUE) %>% 
    ungroup() %>% 
    dplyr::select(-means) %>%
    as.data.frame()
  geneids.to.keep <- dat$gene_id
  rownames(dat) <- dat$gene_symbol
  dat$gene_symbol <- NULL
  dat$gene_id <- NULL
  
  # subset count matrix 
  counts <- counts[which(counts$gene_id %in% geneids.to.keep),]
  counts <- as.data.frame(counts)
  rownames(counts) <- counts$gene_symbol
  counts$gene_id <- NULL
  counts$gene_symbol <- NULL
  
  # do some formatting
  rownames(meta) <- meta$id
  meta <- meta[order(meta$type),]
  dat <- dat[,rownames(meta)]
  counts <- counts[,rownames(meta)]
  meta$id <- NULL
  meta$type <- as.character(meta$type)
  meta$type[meta$type == "Cells - EBV-transformed lymphocytes"] <- "EBV transformed lymphocytes"
  meta$type <- gsub(' ','_',meta$type)
  
  if(identical(rownames(meta), colnames(dat)) & identical(rownames(meta), colnames(counts))){
    print("Dimensions match")
  }
  
  # save full expression matrix
  save(dat, meta, file = 'data/expression_data/gtex_cdls_matrix.RData')
  save(counts, meta, file = 'data/expression_data/gtex_cdls_matrix_counts.RData')
  
  # save mendelian genes expression matrix
  mend.genes <- read.delim('data/Mendelian_Disorders_genelist.txt', stringsAsFactors = F, header = F)
  mend.genes.mat <- dat[which(rownames(dat) %in% mend.genes$V1),]
  mend.genes.counts <- counts[which(rownames(counts) %in% mend.genes$V1),]
  save(mend.genes.mat, mend.genes.counts, meta, file = 'data/expression_data/mendelian_genes_matrix.RData')
}

# pheatmap
meta$id <- rownames(meta) 
meta <- meta[order(meta$type),]
mend.genes.mat <- mend.genes.mat[,rownames(meta)]
meta$id <- NULL
meta$type <- as.character(meta$type)

if(identical(rownames(meta), colnames(mend.genes.mat))){
  print("Dimensions match")
}

# convert to z-scores
myZ <- function(x){
  (x-mean(x))/sd(x)
}
mend.genes.mat.zscored <- data.frame(t(apply(log2(mend.genes.mat + 1), FUN=myZ, MARGIN=1)), check.names = FALSE)

# this shows that expression is not all that different between the three groups
p <- pheatmap(mat = mend.genes.mat.zscored, main = "NMD Genes z-scores",
              annotation_col = meta, 
              show_rownames = FALSE, cluster_rows = F, show_colnames = FALSE)
ggsave(filename = "paper/mendelian_disorders_genes/nmd_genes_zscores_heatmap.png", height = 6, width = 8, device = "png", plot = p)

# normal heatmap with FPKM values and z-score
p <- pheatmap(mat = mend.genes.mat, main = "NMD Genes FPKM",
              annotation_col = meta, 
              show_rownames = FALSE, 
              cluster_rows = F, show_colnames = FALSE, scale = 'row')
ggsave(filename = "paper/mendelian_disorders_genes/nmd_genes_heatmap.png", height = 6, width = 8, device = "png", plot = p)

# check if any of the mendelian genes are in the differential expression list (they shouldn't be)