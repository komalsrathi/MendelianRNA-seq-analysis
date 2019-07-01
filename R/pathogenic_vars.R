# Author: Komal S Rathi
# Date: 02/14/2019
# Function: script to check how many pathogenic variants are caught
# by either one or a combination of callers

library(dplyr)
setwd('~/Projects/DGD_Mendelian_RNASeq/')
final <- read.delim('data/variant_filtering/final_splicevariants_exons.txt', stringsAsFactors = F)

# gatk4
gatk4 <- read.delim('data/variant_filtering/gatk4_filtered_variants.maf')
gatk4 <- gatk4[which(gatk$Tumor_Sample_Barcode %in% final$Sample & gatk4$HGVSc %in% final$HGVSc),]
gatk4 <- gatk4[,c('Tumor_Sample_Barcode','HGVSc','Variant_Classification')]
gatk4$caller <- 'GATK4'

# gatk3.8
gatk <- read.delim('data/variant_filtering/gatk3_filtered_variants.maf')
gatk <- gatk[which(gatk$Tumor_Sample_Barcode %in% final$Sample & gatk$HGVSc %in% final$HGVSc),]
gatk <-gatk[,c('Tumor_Sample_Barcode','HGVSc','Variant_Classification')]
gatk$caller <- 'GATK3.8'

# vardict
vardict <- read.delim('data/variant_filtering/vardict_filtered_variants.maf')
vardict <- vardict[which(vardict$Tumor_Sample_Barcode %in% final$Sample & vardict$HGVSc %in% final$HGVSc),]
vardict <- vardict[,c('Tumor_Sample_Barcode','HGVSc','Variant_Classification')]
vardict$caller <- 'Vardict'

# strelka
strelka <- read.delim('data/variant_filtering/strelka_filtered_variants.maf')
strelka <- strelka[which(strelka$Tumor_Sample_Barcode %in% final$Sample & strelka$HGVSc %in% final$HGVSc),]
strelka <- strelka[,c('Tumor_Sample_Barcode','HGVSc','Variant_Classification')]
strelka$caller <- 'strelka'

# total
total <- rbind(gatk, gatk4, vardict, strelka)
total <- total %>% group_by(Tumor_Sample_Barcode, HGVSc, Variant_Classification) %>%
  summarise(caller = toString(caller), 
         caller.count = n())
total$Sure <- 'Yes'
total[which(total$Tumor_Sample_Barcode %in% c('27571P','95-0614-P','95-1682-P')),'Sure'] <- 'No'
write.table(total, file = 'paper/variant_plots/pathogenic_variants.txt', quote = F, sep = "\t", row.names = F)
