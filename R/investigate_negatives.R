# Author: Komal S. Rathi
# Date: 08/13/2018
# Function: Investigate negative samples that have splice variants

setwd('~/Projects/DGD_Mendelian_RNASeq/')
library(data.table)

# CDL-679 Exon 33 and Exon 34 skipping
# gatk
one <- fread('results/negatives_slides/CDL-679-14P-gatk-haplotype-annotated.maf', skip = 1)
tmp.one <- one[which(one$Hugo_Symbol == "NIPBL"),]
tmp.one <- tmp.one[which(tmp.one$Start_Position >= 37027268 & tmp.one$End_Position <= 37044650),]
# vardict
two <- fread('results/negatives_slides/Sample_1__CDL-679-14P.vardict-annotated-rnaedit-annotated-gemini.maf', skip = 1)
tmp.two <- two[which(two$Hugo_Symbol == "NIPBL"),]
tmp.two <- tmp.two[which(tmp.two$Start_Position >= 37027268 & tmp.two$End_Position <= 37044650),]
write.table(tmp.one[,c(6:12,14,35:37,40:42)], file = 'results/negatives_slides/CDL-679_potential_variants.txt', quote = F, sep = "\t", row.names = F)

# 95-0614
# gatk
one <- fread('results/negatives_slides/95-0614-P-gatk-haplotype-annotated.maf', skip = 1)
tmp.one <- one[which(one$Hugo_Symbol == "NIPBL"),]
tmp.one <- tmp.one[which(tmp.one$Start_Position >= 37010152 & tmp.one$End_Position <= 37014836),]
# vardict
two <- fread('results/negatives_slides/Sample_1__95-0614-P.vardict-annotated-rnaedit-annotated-gemini.maf', skip = 1)
tmp.two <- two[which(two$Hugo_Symbol == "NIPBL"),]
tmp.two <- tmp.two[which(tmp.two$Start_Position >= 37010152 & tmp.two$End_Position <= 37014836),]
write.table(tmp.one[,c(6:12,14,35:37,40:42)], file = 'results/negatives_slides/95-0614-P_potential_variants.txt', quote = F, sep = "\t", row.names = F)
