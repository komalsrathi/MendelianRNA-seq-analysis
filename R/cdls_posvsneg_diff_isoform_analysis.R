# Author: Komal S. Rathi
# Date: 05/02/2019
# Function: Differential Isoform Analysis 
# between negative (CDL-022, CDL-069, CDL-086) and positive samples

setwd('~/Projects/DGD_Mendelian_RNASeq/')
library(reshape2)
library(limma)
library(ggrepel)
source('R/pubTheme.R')

negative <- 'CDL-219-P'
dat <- readRDS('data/expression_data/CDLS_patient_expr_transcripts_hg38.RDS')
dat <- dat[which(dat$Sample != negative),]
samples <- unique(dat$Sample)

# cdls genes
cdls <- read.delim('data/CDL_genes.txt', header = F, stringsAsFactors = F)
cdls <- cdls$V1

# expr matrix
expr <- dcast(dat, transcript_id~Sample, value.var = 'expected_count')
rownames(expr) <- expr$transcript_id
expr$transcript_id <- NULL
expr.sub <- expr[apply(expr, MARGIN = 1, FUN = function(x) all(x >= 10)),]
expr.sub <- expr.sub[,order(colnames(expr.sub))]

# metadata
metadata <- data.frame(Samples = samples, type = '')
metadata$type <- ifelse(metadata$Samples %in% c("CDL-022-P", "CDL-069-99P", "CDL-086-99P"), "negative", "positive")
rownames(metadata) <- metadata$Samples
metadata <- metadata[order(rownames(metadata)),]

# pca (no clustering of groups observed)
pcadat <- prcomp(log2(expr.sub+1), scale. = T, center = T)
pcadat <- pcadat$rotation
pcadat <- data.frame(pcadat)[1:4]
pcadat <- merge(pcadat, metadata, by.x = 'row.names', by.y = 'Samples')
p <- ggplot(pcadat, aes(x = PC1, y = PC2, label = Row.names)) + 
  # geom_label(aes(fill = factor(type)), colour = "white", fontface = "bold") +
  geom_point(aes(fill = type), pch = 21, size = 3) +
  geom_text_repel() +
  theme_bw() +
  theme_Publication2() + 
  theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  ggtitle(paste0("PCA\nExpressed Transcripts\n(n = ",nrow(expr.sub),")"))
p

# limma analysis
if(identical(rownames(metadata), colnames(expr.sub))){
  print("ok")
}
mydesign <- model.matrix(~0+type, metadata)
colnames(mydesign) <- gsub('type','',colnames(mydesign))
if(identical(rownames(mydesign), colnames(expr.sub))){
  print("ok")
}
tmp.voom <- voom(counts = as.matrix(expr.sub), design = mydesign)
tmp.voom <- tmp.voom$E
fit <- lmFit(tmp.voom, design = mydesign)
contrast.matrix = makeContrasts(contrasts = "positive-negative", levels = mydesign)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res <- topTable(fit2, number=Inf, 1)[,c("logFC", "P.Value", "adj.P.Val")]

# merge to gene annotation
dat <- data.table::fread('~/Projects/toil-rnaseq-20k/data/annotation/gencode.v23.annotation.txt')
dat <- unique(dat[,c('gene_symbol','transcript_id')])
dat <- as.data.frame(dat)
dat <- merge(dat, res, by.x = 'transcript_id', by.y = 'row.names')
dat <- dat[which(dat$gene_symbol %in% cdls),]

# nothing came out to be significant because individual groups do not cluster within themselves