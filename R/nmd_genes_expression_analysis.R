# Author: Komal S. Rathi
# Date: 04/15/2019
# Function: 
# Filter out low expressing mendelian genes (FPKM < 1 and count < 10)
# Barplot and Boxplot
# Do Venn
# Do PCA

setwd('~/Projects/DGD_Mendelian_RNASeq/')
library(ggpubr)
library(reshape2)
library(dplyr)
library(scatterplot3d)
library(VennDiagram)
library(stringi)
require(gridExtra)

source('R/pubTheme.R')

# annotation - only protein coding genes
# total 19797 protein coding genes
annot <- data.table::fread('~/Projects/toil-rnaseq-20k/data/annotation/gencode.v23.annotation.txt')
annot <- unique(annot[which(annot$biotype == "protein_coding"),c('gene_id','gene_symbol')])
annot <- as.data.frame(annot)
annot$ensembl <- gsub('[.].*','',annot$gene_id)

# patient raw counts
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

if(file.exists('data/mend_genes_isoform_analysis_input.RData')){
  print('Exists')
  load('data/mend_genes_isoform_analysis_input.RData')
} else {
  p.expr <- expr.genes(rds.path = 'data/expression_data/CDLS_patient_expr_hg38.RDS')
  wb.expr <- expr.genes(rds.path = 'data/expression_data/WB_expr_hg38.RDS')
  ebv.expr <- expr.genes(rds.path = 'data/expression_data/EBV_expr_hg38.RDS')
  
  # full expression matrix
  load('data/expression_data/gtex_cdls_matrix.RData')
  meta$id <- rownames(meta)
  
  patient <- dat[,colnames(dat) %in% meta[which(meta$type == "CdLS_Patient"),'id']]
  patient <- melt(as.matrix(patient), varnames = c("gene_symbol","id"), value.name = "FPKM")
  patient <- patient[which(patient$gene_symbol %in% annot$gene_symbol),]
  
  # filter out low expressing protein coding genes from patient cohort
  # we use average FPKM > 1 as cutoff
  # 11243 genes are expressed with average FPKM > 1
  patient <- patient %>% group_by(gene_symbol) %>%
    summarise(FPKM = mean(FPKM)) %>%
    filter(FPKM > 1) %>%
    as.data.frame()
  
  # how many of these are mendelian genes
  mend.genes <- read.delim('data/Mendelian_Disorders_genelist.txt', stringsAsFactors = F, header = F)
  patient$mend.gene <- ifelse(patient$gene_symbol %in% mend.genes$V1, 'Y', 'N')
  patient <- patient[which(patient$gene_symbol %in% p.expr$gene_symbol),]
  plyr::count(patient$mend.gene) # 9299 (non-mendelian), 1745 (mendelian), 11044 all expressed
  patient$type <- 'CdLS_Patients'
  
  # Out of the 652 genes that are not expressed, let's see how many are expressed in GTEx tissues:
  gtex <- dat[,colnames(dat) %in% meta[which(meta$type != "CdLS_Patient"),'id']]
  gtex <- melt(as.matrix(gtex), varnames = c("gene_symbol","id"), value.name = "FPKM")
  gtex <- gtex[which(gtex$gene_symbol %in% annot$gene_symbol),]
  gtex <- merge(gtex, meta, by = 'id')
  gtex <- gtex %>% group_by(gene_symbol, type) %>%
    summarise(FPKM = mean(FPKM)) %>%
    filter(FPKM > 1) %>%
    as.data.frame()
  gtex$mend.gene <- ifelse(gtex$gene_symbol %in% mend.genes$V1, 'Y', 'N')
  gtex.wb <- gtex[which(gtex$type == "Whole_Blood"),] # expressed in WB
  gtex.wb <- gtex.wb[which(gtex.wb$gene_symbol %in% wb.expr$gene_symbol),]
  gtex.ebv <- gtex[which(gtex$type == "EBV_transformed_lymphocytes"),] # expressed in EBV
  gtex.ebv <- gtex.ebv[which(gtex.ebv$gene_symbol %in% ebv.expr$gene_symbol),]
  
  # save these big inputs
  save(gtex.ebv, gtex.wb, patient, p.expr, ebv.expr, wb.expr, file = 'data/mend_genes_isoform_analysis_input.RData')
}

# barplot of Mendelian Neuro genes with FPKM > 1 in EBV vs Whole blood
for.barplot <- rbind(gtex.ebv, gtex.wb, patient) 
for.barplot <- for.barplot[which(for.barplot$mend.gene == "Y"),]
for.boxplot <- dcast(for.barplot, gene_symbol~type, value.var = 'FPKM', fill = 0)
for.boxplot <- melt(for.boxplot, value.name = 'FPKM', variable.name = 'Groups')
my_comparisons <- list(c("EBV_transformed_lymphocytes", "Whole_Blood"),
                       c("CdLS_Patients", "Whole_Blood"),
                       c("CdLS_Patients", "EBV_transformed_lymphocytes"))
p <- ggplot(for.boxplot, aes(x = Groups, y = log2(FPKM + 1), fill = Groups)) + 
  stat_boxplot(geom ='errorbar', width = 0.2) +
  geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) + 
  theme_Publication() +
  ggtitle("NMD Genes Average Expression\nFPKM > 1 & Coverage > 10") + guides(fill = FALSE) +
  stat_compare_means(comparisons = my_comparisons) + ylab('log2(FPKM)')
ggsave(filename = 'paper/mendelian_disorders_genes/expressed_nmd_genes_boxplot.png', plot = p, device = 'png', width = 7, height = 6)

for.barplot <- plyr::count(for.barplot, c("type","mend.gene"))
for.barplot$total <- 2541
for.barplot$prop <- round(for.barplot$freq/for.barplot$total*100, 2)
for.barplot$label <- paste0(for.barplot$freq,' (',for.barplot$prop,'%)')
p <- ggplot(for.barplot, aes(x = type, y = prop, fill = type)) + 
  geom_bar(stat = 'identity', color = 'black',  width=0.8) +
  theme_bw() + 
  theme_Publication2() + ylab('% of Expressed Mendelian Genes') + xlab('') +
  ggtitle('Proportion of Expressed NMD Genes (%)\nFPKM > 1 & Coverage > 10') + 
  geom_text(aes(label = label, vjust = 1.5, hjust = 0.5), position = position_dodge(width = 0.8), size = 3) +
  guides(fill = FALSE)
p
ggsave(filename = 'paper/mendelian_disorders_genes/expressed_nmd_genes_barplot.png', plot = p, device = 'png', width = 6, height = 4)

# more analysis (revise this part)
# gtex <- dcast(gtex, gene_symbol + mend.gene ~ type, value.var = 'FPKM')

# EBV
# plyr::count(gtex.ebv$mend.gene) # 8906 (non-mendelian), 1706 (mendelian), 10612 all expressed

# WB
# plyr::count(gtex.wb$mend.gene) # 4700 (non-mendelian), 917 (mendelian), 5617 all expressed 

# venn diagram of expressed protein coding genes
png(filename = 'paper/expression_analysis/expressedgenes_overlap_venn.png', width = 6, height = 5, units = "in", res = 240)
one <- patient$gene_symbol
two <- gtex.ebv$gene_symbol
three <- gtex.wb$gene_symbol
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
                                   paste0("WB\n(",length(three),")")),
                      fill = c("skyblue", "pink1", "mediumorchid"), lty = "blank")
grid.arrange(gTree(children=p), top="Expressed genes")
dev.off()
# ggsave(filename = 'paper/expression_analysis/expressedgenes_overlap_venn.png', device = 'png', plot = p, width = 7, height = 5)

# venn diagram of expressed mendelian genes
png(filename = 'paper/mendelian_disorders_genes/expressed_nmd_genes_overlap_venn.png', width = 6, height = 5, units = "in", res = 240)
one <- patient[which(patient$mend.gene == 'Y'),'gene_symbol']
two <- gtex.ebv[which(gtex.ebv$mend.gene == 'Y'),'gene_symbol']
three <- gtex.wb[which(gtex.wb$mend.gene == 'Y'),'gene_symbol']
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
                                   paste0("WB\n(",length(three),")")),
                      fill = c("skyblue", "pink1", "mediumorchid"), lty = "blank")
grid.arrange(gTree(children=p), top="Expressed NMD genes")
dev.off()
# ggsave(filename = 'paper/mendelian_disorders_genes/mendeliangenes_overlap_venn.png', device = 'png', plot = p, width = 7, height = 5)

# do a PCA of expressed genes (in all three types)
one <- patient$gene_symbol
two <- gtex.ebv$gene_symbol
three <- gtex.wb$gene_symbol
expr.genes <- intersect(intersect(one, two), three)
dat <- dat[which(rownames(dat) %in% expr.genes),]

pcadat <- prcomp(log2(dat+1), scale. = T, center = T)
pcadat <- pcadat$rotation
pcadat <- data.frame(pcadat)[1:4]
pcadat <- merge(pcadat, meta, by.x = 'row.names', by.y = 'id')
pcadat$color <- ifelse(pcadat$type == "CdLS_Patient", "#00AFBB", ifelse(pcadat$type == "Whole_Blood", "#FC4E07", "#E7B800"))
pcadat$shape <- ifelse(pcadat$type == "CdLS_Patient", 0, ifelse(pcadat$type == "Whole_Blood", 1, 2))

# 2D
p <- ggplot(pcadat, aes(x = PC1, y = PC2, colour = type)) + 
  geom_point() + 
  theme_bw() +
  theme_Publication2() + 
  theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  scale_color_manual(name = 'Type', values = c("CdLS_Patient" = "#00AFBB", 
                                "Whole_Blood" = "#FC4E07",
                                "EBV_transformed_lymphocytes" = "#E7B800")) +
  ggtitle(paste0("PCA\nExpressed Protein coding genes\n(n = ",nrow(dat),")"))
p
ggsave(plot = p, filename = 'paper/expression_analysis/PCA_expressed_genes_2D.png', device = 'png', height = 6, width = 8)

# 3D
fname <- paste0('paper/expression_analysis/PCA_expressed_genes_3D.png')
png(filename = fname, width = 1000, height = 800, res = 120)
s3d <- scatterplot3d(pcadat[,"PC2"], pcadat[,"PC3"], pcadat[,"PC4"], 
              xlab="PC2", ylab="PC3", zlab="PC4", 
              color=as.character(pcadat[,"color"]),
              main=paste0("PCA\nExpressed Protein coding genes\n(n = ",nrow(dat),")"), pch = 20,
              cex.symbols=1)
legend(x = 'bottomright', legend = unique(pcadat$type), 
       pch = 16, col =  c("#00AFBB", "#FC4E07", "#E7B800"), 
       cex = 0.7, inset=0)
dev.off()

# genes that have at least 1 patient with FPKM > 1 
# mendgenes.expr <- patient[which(patient$mend.gene == "Y"),]
# not.expressed <- setdiff(mend.genes$V1, mendgenes.expr$gene_symbol)
# write.table(not.expressed, file = 'paper/mendelian_disorders_genes/mend_genes_notexpressed_cdls.txt', quote = F, sep = "\t", row.names = F, col.names = F)
# write.table(mendgenes.expr$gene_symbol, file = 'paper/mendelian_disorders_genes/mend_genes_expressed_cdls.txt', quote = F, sep = "\t", row.names = F, col.names = F)

# do GO/pathway analysis on expressed mendelian genes
