# Author: Komal S. Rathi
# Date: 06/14/2019
# Function: 
# Differential Expression Analysis between all three groups
# Identify the proportion of Differentially Expressed Mendelian Genes

library(limma)
library(reshape2)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(ggpubr)

setwd('~/Projects/DGD_Mendelian_RNASeq')
source('R/pubTheme.R')

if(file.exists('data/expression_data/diffexpr_limma.RData')){
  print("Exists")
  load('data/expression_data/diffexpr_limma.RData')
} else {
  # dataset (FPKM values)
  load('data/expression_data/gtex_cdls_matrix.RData')
  # dataset (counts)
  load('data/expression_data/gtex_cdls_matrix_counts.RData')
  dat <- counts
  
  if(identical(rownames(meta), colnames(dat))){
    print("Dimensions match")
  }
  
  # remove all zeros
  dat <- dat[rowSums(dat) > 0,]
  
  # normalize the dataset
  mydesign <- model.matrix(~0+type, meta)
  colnames(mydesign) <- gsub('type','',colnames(mydesign))
  dat.voom <- voom(counts = as.matrix(dat), design = mydesign)
  dat.voom <- dat.voom$E
  
  # do differential expression between all three groups
  fit <- lmFit(dat.voom, design = mydesign)
  contrast.matrix = makeContrasts(contrasts = c('CdLS_Patient-EBV_transformed_lymphocytes',
                                                'CdLS_Patient-Whole_Blood',
                                                'EBV_transformed_lymphocytes-Whole_Blood'),
                                  levels = mydesign)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  cdls.vs.ebv <- topTable(fit2, number=Inf, 1)
  cdls.vs.wb <- topTable(fit2, number = Inf, 2)
  ebv.vs.wb <- topTable(fit2, number = Inf, 3)
  
  # what proportion of differentially expressed genes overlap between contrast 2 and 3.
  # scatterplot of overlapping genes 
  # ebv.vs.wb and cdls.vs.wb
  add.annot <- function(x, txt){
    x$DEGene <- x$adj.P.Val < 0.05 
    x$Label <- ifelse(x$DEGene == TRUE & x$logFC < 0,'Downreg', 'Upreg')
    x$Label <- ifelse(x$DEGene == FALSE,'NotSig',x$Label)
    x$Comparison <- txt
    colnames(x) <- paste0(txt, '_', colnames(x))
    return(x)
  }
  cdls.vs.ebv <- add.annot(cdls.vs.ebv,'CdLS_vs_EBV')
  cdls.vs.wb <- add.annot(cdls.vs.wb,'CdLS_vs_WB')
  ebv.vs.wb <- add.annot(ebv.vs.wb,'EBV_vs_WB')
  
  # save gene lists
  save(cdls.vs.ebv, cdls.vs.wb, ebv.vs.wb, file = 'data/expression_data/diffexpr_limma.RData')
}

total <- merge(cdls.vs.wb[,c(1,5,8)], ebv.vs.wb[,c(1,5,8)], by = 'row.names')
total$NewLabel <- NA
total$NewLabel[total$CdLS_vs_WB_Label == "Upreg" & total$EBV_vs_WB_Label == "Upreg"] <- 'Upregulated'
total$NewLabel[total$CdLS_vs_WB_Label == "Downreg" & total$EBV_vs_WB_Label == "Downreg"] <- 'Downregulated'
total$NewLabel[is.na(total$NewLabel)] <- 'Others'
nums <- plyr::count(total$NewLabel)
total <- merge(total, nums, by.x = 'NewLabel', by.y = 'x')
total$NewLabel <- paste0(total$NewLabel," (n = ", total$freq, ")")

summ <- summary(lm(formula = CdLS_vs_WB_logFC ~ EBV_vs_WB_logFC, data = total))
summ

# scatter plot
ggplot(total, aes(CdLS_vs_WB_logFC, EBV_vs_WB_logFC)) +
  geom_point(size = 0.5, aes(color = NewLabel)) +
  xlab('CdLS vs WB (logFC)') + 
  ylab('EBV vs WB (logFC)') + theme_bw() +
  theme_Publication() + 
  ggtitle('Correlation of logFC\n CdLS vs Whole Blood and EBV-Lymphocytes vs Whole Blood') +
  guides(color=guide_legend(title="Expression Status")) + 
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
  annotate("text", x=5, y=-5, label = "R^2 = 0.8909\nP-value < 2.2e16") + 
  scale_color_manual(values=c("red", "gray", "#00BA38"))
ggsave("paper/diffexpr_comparisons/CdLSAndEBV_vs_WB_logFC_comparison.png", device = 'png', width=10, height=8)

# modify data for dodged barplot of mendelian genes
forbar <- merge(cdls.vs.ebv[,c(8,9)], cdls.vs.wb[,c(8,9)], by = 'row.names')
forbar <- merge(forbar, ebv.vs.wb[,c(8,9)], by.x = 'Row.names', by.y = 'row.names')
forbar <- forbar[,grep('Row.names|Label', colnames(forbar))]
colnames(forbar) <- gsub('_Label', '', colnames(forbar))
forbar <- melt(forbar, id.vars = 'Row.names')

# add mendelian genes
mend.genes <- read.delim('data/Mendelian_Disorders_genelist.txt', stringsAsFactors = F, header = F)
mend.genes <- mend.genes$V1
forbar$Mendelian_Genes[forbar$Row.names %in% mend.genes] <- 'Yes'
forbar$Mendelian_Genes[is.na(forbar$Mendelian_Genes)] <- 'No'
ct <- plyr::count(forbar, c("variable","Mendelian_Genes","value"))
ct <- ct %>% group_by(variable, value) %>% 
  mutate(total = sum(freq)) %>%
  mutate(perc = round(freq/total*100, 2))

# dodged barplot of mendelian genes
ct$Mendelian_Genes <- factor(ct$Mendelian_Genes, levels = c("Yes","No"))
ct$value <- factor(ct$value, levels = c("Upreg","Downreg","NotSig"))
p <- ggplot(ct, aes(x = value, y = perc, fill = factor(Mendelian_Genes), label = paste0(perc,"%\n(", freq, ")"))) + 
  geom_bar(stat = 'identity', position = 'dodge') + facet_wrap(~variable) +
  labs(x = "Category", y = "Percent (%)", fill = "Mendelian Genes") +
  geom_text(size = 3,
            position = position_dodge(0.9),
            vjust = -0.1) +
  theme_bw() + theme_Publication() + 
  ggtitle('Proportion of Differentially Expressed\nNMD Genes')
ggsave(filename = "paper/diffexpr_comparisons/Proportion_diffexpr_nmd_genes.png", plot = p, device = 'png', width=12, height=8)

# differentially expressed genes
cdls.vs.ebv <- cdls.vs.ebv[which(cdls.vs.ebv$CdLS_vs_EBV_DEGene == TRUE),]
cdls.vs.wb <- cdls.vs.wb[which(cdls.vs.wb$CdLS_vs_WB_DEGene == TRUE),]
ebv.vs.wb <- ebv.vs.wb[which(ebv.vs.wb$EBV_vs_WB_DEGene == TRUE),]

# venn diagram
png(filename = 'paper/diffexpr_comparisons/diffexpr_venn.png', width = 6, height = 5, units = "in", res = 240)
one <- rownames(cdls.vs.ebv)
two <- rownames(cdls.vs.wb)
three <- rownames(ebv.vs.wb)
p <- draw.triple.venn(area1 = length(one), 
                      area2 = length(two), 
                      area3 = length(three), 
                      n12 = length(intersect(one, two)),
                      n23 = length(intersect(two, three)), 
                      n13 = length(intersect(one, three)),
                      n123 = length(intersect(intersect(one, two), three)), 
                      scaled = F, euler.d = F,
                      category = c("CdLS vs EBV", "CdLS vs WB", "EBV vs WB"),
                      fill = c("skyblue", "pink1", "mediumorchid"), lty = "blank")
grid.arrange(gTree(children=p), top="Differentially Expressed Genes")
dev.off()
# ggsave(filename = 'paper/diffexpr_comparisons/venn.png', device = 'png', plot = p, width = 7, height = 5)

# for metacore/pathway analysis
cdls.vs.wb$gene_symbol <- rownames(cdls.vs.wb)
cdls.vs.wb <- cdls.vs.wb[which(cdls.vs.wb$CdLS_vs_WB_DEGene == T & 
                                 cdls.vs.wb$CdLS_vs_WB_adj.P.Val < 0.01 & 
                                 abs(cdls.vs.wb$CdLS_vs_WB_logFC) > 1), c('gene_symbol','CdLS_vs_WB_logFC','CdLS_vs_WB_adj.P.Val')]
write.table(cdls.vs.wb, file = 'data/expression_data/diffexpr_cdls_vs_wb.txt', quote = F, sep = "\t", row.names = F)

