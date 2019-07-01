# Author: Beryl Cummings
# Modified By: Komal S. Rathi
# Date: 04/04/2019
# Function: Determine sex of patients using RNA-seq expression

library(ggplot2)
library(plyr)
library(dplyr)
source('R/sex_qc_part1.R')

gtex_sex_biased_rpkm <- read.delim('data/gtex_expression_sex_biased.txt', stringsAsFactors = F, row.names = 1)
patients_rpkm <- read.delim('data/patient_genes.fpkm.gct', stringsAsFactors = F, skip = 2, header = T, check.names = F)

patients_rpkm <- extractGenesOfInterest(patients_rpkm, sex_biased_genes$gene_id)

genes_in_both <- intersect(names(gtex_sex_biased_rpkm), names(patients_rpkm))
all_sex_biased <- rbind(subset(gtex_sex_biased_rpkm, select = genes_in_both), subset(patients_rpkm, select = genes_in_both))

PCADat <- performPCA(all_sex_biased, "sex",expected_annotations =c("male","female"))
PCADat$sex[PCADat$sex == "female"] <- "female_GTEx"
PCADat$sex[PCADat$sex == "male"] <- "male_GTEx"

# add known info (Komal)
g.info <- read.delim('data/gender_info.txt', stringsAsFactors = F)
f.patients <- g.info[which(g.info$Gender == 'Female'),1]
PCADat[which(rownames(PCADat) %in% f.patients),'sex'] <- 'female_patient'
PCADat$sex[PCADat$sex == "Patient"] <- "male_patient"

chrom_y_genes <-  subset(sex_biased_genes, chr == "chrY")$gene_id
all_xist_y <- getXISTAvgYChromExpr(all_sex_biased, chrom_y_genes)
all_xist_y[which(all_xist_y$sample %in% f.patients),'sex'] <- 'female_patient'
all_xist_y$sex[all_xist_y$sex == "Patient"] <- "male_patient"
all_xist_y$sex[all_xist_y$sex == "female"] <- "female_GTEx"
all_xist_y$sex[all_xist_y$sex == "male"] <- "male_GTEx"

plot_theme <- theme_bw()+theme(plot.title=element_text(size=20), axis.text=element_text(size=15), axis.title=element_text(size=14), legend.title=element_blank(), legend.text=element_text(size=10),panel.grid.major=element_blank())
colours <- scale_colour_manual(values=c("#700D4F","#EAA8D4","#2f6ba8", "#87B2DD"))


pdfFile = 'paper/expression_analysis/gender_check.pdf'
pdf(pdfFile, width =8, height =6)

# plot(PCA,type="l",main="Variance Explained by PCs")

ggplot( PCADat %>% arrange(sex), aes( x = PC1, y = PC2, col = sex)) + 
  geom_point( size=3  ) + 
  plot_theme + colours +
  xlab( paste( "PC1 (", signif( PoV[1]*100, 2 ), "%)", sep = "" )) + 
  ylab( paste (" PC2 (",signif( PoV[2]*100, 2 ), "%)", sep = ""))

ggplot(all_xist_y %>% arrange(sex), aes( x = XIST, y = avg_y_chrom, col = as.factor(sex))) + 
  geom_point(size = 3)+
  plot_theme + colours + 
  xlab("XIST expression") + ylab("Average of expression of sex-biased genes on y chromosome")

dev.off()

# hclustfunc <- function(x) hclust(x, method="complete")
# distfunc <- function(x) as.dist(1-cor(t(x),method = "spearman"))
# d <- distfunc(all_sex_biased)
# fit <- hclustfunc(d)
# allClusters<-cutree(fit, k = 2)

# for (clusnum in 1:2){
#   clus <-  allClusters[allClusters==clusnum]
#   cluster <- unique(gsub("\\.[0-9]*$","",names(clus)))
#   cat('\n')
#   if('male' %in% cluster){
#     cat("\nMale patients are:")
#     for(elem in names(clus)){
#       if(grepl("male",elem)){next}
#       cat(paste("\n", elem))
#     }
#   }
#   if('female' %in% cluster){
#     cat("\nFemale patients are:")
#     for(elem in names(clus)){
#       if(grepl("female",elem)){next}
#       cat(paste("\n", elem))
#     }
#   }
# }
