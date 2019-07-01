# Author: Komal S. Rathi
# Date: 12/23/2018
# Function: Average splice events per sample summary in Patient and GTex samples

setwd('~/Projects/DGD_Mendelian_RNASeq/')
library(reshape2)
library(dplyr)
library(tidyr)

patient <- read.delim('results/splice-events-patient/All.CDL_genes.normalized.splicing.txt', stringsAsFactors = F)
gtex <- read.delim('results/splice-events-gtex/All.CDL_genes.normalized.splicing.txt', stringsAsFactors = F)

patient.full <- patient %>% 
  mutate(Samples = strsplit(as.character(Samples), ",")) %>% 
  unnest(Samples)
patient.full$Samples <- gsub('.*:','',patient.full$Samples)
patient.full <- unique(patient.full[,c('Gene','Pos','Samples')])
patient.ct <- plyr::count(patient.full, c("Samples","Gene"))
patient.ct <- patient.ct %>% group_by(Samples) %>% 
  summarise(ave.variants.per.gene = round(mean(freq),2), total.variants.per.sample = sum(freq))
write.table(patient.ct, file = 'results/splice-events-summary/patient_splice_events.txt', quote = F, sep = "\t", row.names = F)

gtex.full <- gtex %>%
  mutate(Samples = strsplit(as.character(Samples), ",")) %>% 
  unnest(Samples)
gtex.full$Samples <- gsub('.*:','',gtex.full$Samples)
gtex.full <- unique(gtex.full[,c('Gene','Pos','Samples')])
gtex.ct <- plyr::count(gtex.full, c("Samples","Gene"))
gtex.ct <- gtex.ct %>% group_by(Samples) %>% 
  summarise(ave.variants.per.gene = round(mean(freq),2), total.variants.per.sample = sum(freq))
write.table(gtex.ct, file = 'results/splice-events-summary/normal_splice_events.txt', quote = F, sep = "\t", row.names = F)
