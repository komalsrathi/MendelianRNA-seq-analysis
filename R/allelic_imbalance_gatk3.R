# Author: Komal S. Rathi
# Date: 04/19/2018
# Function: Compute Allelic Imbalance (specifically for HDAC8 in CDL-418)

setwd('/mnt/isilon/cbmi/variome/rathik/mendelian_rnaseq/')

library(data.table)
library(dplyr)

# gatk3.8
lf <- list.files(path = 'updated-mafs/gatk3-source-maf/', pattern = "hgmdannotated.maf", full.names = TRUE)
for(i in 1:length(lf)){
  print(lf[i])
  dat <- fread(lf[i])
  dat$Tumor_Sample_Barcode <- gsub('.*/|-hgmdannotated.maf', '', lf[i])
  dat <- dat[which(dat$t_depth >= 20),] # filter by depth >= 20 (change to 10)
  dat$AB_ref <- dat$t_ref_count/dat$t_depth # allele balance (ref)
  dat$AB_alt <- dat$t_alt_count/dat$t_depth # allele balance (alt)
  dat <- dat %>% group_by(Hugo_Symbol) %>% # calculate number of variants per gene
    mutate(var.count = n()) %>%
    as.data.frame()
  dat <- dat %>% filter(VARIANT_CLASS == "SNV", 
                        HGVSc != "") %>% # only keep heterozygous SNV calls
    as.data.frame()
  dat <- unique(dat[,c("Hugo_Symbol","Reference_Allele",
    "Tumor_Seq_Allele1","Tumor_Seq_Allele2","t_ref_count","t_alt_count","t_depth",
    "HGVSc","ALLELE_NUM","VARIANT_CLASS",
    "AB_ref","AB_alt","var.count","Tumor_Sample_Barcode")]) 
  if(i == 1){
    total <- dat
  } else {
    total <- rbind(total, dat)
  }
}
gatk3 <- total
save(gatk3, file = 'allelic-imbalance/gatk3_allelic_imbalance.RData')
