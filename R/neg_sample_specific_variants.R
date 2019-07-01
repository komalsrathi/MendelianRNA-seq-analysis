# Author: Komal S Rathi
# Date: 05/20/2019
# Function: Script to keep unique variants called in negative samples

setwd('~/Projects/DGD_Mendelian_RNASeq/')

library(dplyr)

# all samples gatk total
gatk.lf <- list.files('data/variant_filtering/rawdata/gatk4', pattern = '*hgmdannotated.maf', full.names = T)
gatk.lf <- grep('CDL-219', gatk.lf, invert = T, value = T)
for(i in 1:length(gatk.lf)){
  print(i)
  if(i == 1){
    gatk.total <- data.table::fread(gatk.lf[i])
    gatk.total$Tumor_Sample_Barcode <- gsub('.*/|-gatk.*','',gatk.lf[i])
  } else {
    dat <- data.table::fread(gatk.lf[i])
    dat$Tumor_Sample_Barcode <- gsub('.*/|-gatk.*','',gatk.lf[i])
    gatk.total <- rbind(gatk.total, dat)
  }
}
gatk.total <- unique(gatk.total[,c('Hugo_Symbol','HGVSc','Start_Position','End_Position','Tumor_Sample_Barcode')])
gatk.total <- gatk.total %>% 
  group_by(Hugo_Symbol,HGVSc,Start_Position,End_Position) %>%
  summarise(x = n()) %>% filter(x == 15) %>%
  unique() %>%
  as.data.frame()

# all samples vardict total
vardict.lf <- list.files('data/variant_filtering/rawdata/vardict', pattern = '*hgmdannotated.maf', full.names = T)
vardict.lf <- grep('CDL-219', vardict.lf, invert = T, value = T)

for(i in 1:length(vardict.lf)){
  print(i)
  if(i == 1){
    vardict.total <- data.table::fread(vardict.lf[i])
    vardict.total$Tumor_Sample_Barcode <- gsub('.*/|.vardict.*','',vardict.lf[i])
  } else {
    dat <- data.table::fread(vardict.lf[i])
    dat$Tumor_Sample_Barcode <- gsub('.*/|.vardict.*','',vardict.lf[i])
    vardict.total <- rbind(vardict.total, dat)
  }
}
vardict.total <- unique(vardict.total[,c('Hugo_Symbol','HGVSc','Start_Position','End_Position','Tumor_Sample_Barcode')])
vardict.total <- vardict.total %>% 
  group_by(Hugo_Symbol,HGVSc,Start_Position,End_Position) %>%
  summarise(x = n()) %>% filter(x == 15) %>%
  unique() %>%
  as.data.frame()

to.remove <- unique(rbind(vardict.total, gatk.total))

# mutations from three 3 negatives
neg <- c("CDL-022-P", "CDL-069-99P", "CDL-086-99P")
paste0(neg,collapse = "|")

# GATK result (F3) for negatives
gatk <- read.delim('data/variant_filtering/gatk4_filtered_variants.maf')
gatk <- gatk[grep(paste0(neg,collapse = "|"), gatk$Tumor_Sample_Barcode),]
gatk$caller <- 'gatk'

# Vardict result (F4) for negatives
vardict <- read.delim('data/variant_filtering/vardict_filtered_variants.maf')
vardict <- vardict[grep(paste0(neg,collapse = "|"), vardict$Tumor_Sample_Barcode),]
vardict$caller <- 'vardict'

# total
total <- rbind(gatk, vardict)
total <- total[-which(total$Hugo_Symbol %in% to.remove$Hugo_Symbol &
              total$HGVSc %in% to.remove$HGVSc &
              total$Start_Position %in% to.remove$Start_Position &
              total$End_Position %in% to.remove$End_Position),]
total <- total[order(total$Tumor_Sample_Barcode, total$Hugo_Symbol, total$caller),]
total <- total[,colSums(is.na(total))<nrow(total)]
write.table(total, file = 'paper/variant_plots/variants_unique_neg_samples.txt', quote = F, sep = "\t", row.names = F)
