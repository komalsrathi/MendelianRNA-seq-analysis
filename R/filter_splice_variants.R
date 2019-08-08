# Author: Komal S. Rathi
# Date: 12/23/2018
# Function: Filter and summarize disruptive splice variants
# inspired from Cummings et al.

setwd('~/Projects/DGD_Mendelian_RNASeq/')
library(reshape2)

# read data
patient <- read.delim('results/splice-events-patient/All.CDL_genes.normalized.splicing.txt', stringsAsFactors = F)
gtex <- read.delim('results/splice-events-gtex/All.CDL_genes.normalized.splicing.txt', stringsAsFactors = F)
total <- data.frame()

############################################
# n = 6 (all six events)
# Method 1: Events that are only seen in one patient sample, with high read support.
# So, here all the events that were found only in 1 patient with high read support were recovered because they have very few reads (n=1) supporting in gtex samples.
only1sample <- read.delim('results/splice-events-patient/FilteredJunctions_ReadSupport10_Only1Sample.txt')
only1sample.gtex <- gtex[which(gtex$Pos %in% only1sample$Pos),] # all low-events so recover everything from only1sample
total <- rbind(total,only1sample)
############################################

############################################
# Method 2: Events seen in many individuals, but only seen in one individual with read support higher than say, 10 reads (we can make it more stringent i.e. 20 reads)
# filter < 5 reads (n = 9)
filtered.reads.5 <- read.delim('results/splice-events-patient/FilteredJunctions_ReadSupport10_FilterReads5.txt')
filtered.reads.5 <- filtered.reads.5[which(filtered.reads.5$NSamplesSeen == 1),]
filtered.reads.5 <- filtered.reads.5[which(filtered.reads.5$NTimesSeen >= 10),]
filtered.reads.5 <- filtered.reads.5[-which(filtered.reads.5$Pos %in% only1sample$Pos),]

# n = 3
total <- rbind(total, filtered.reads.5)
total$PropSeen.Samples <- gsub(':.*','',total$PropSeen.Samples)
total$PropSeen.Samples <- as.numeric(total$PropSeen.Samples)*100
total <- unique(total[,c('Sample','Gene','Pos','Annotation','NTimesSeen','NSamplesSeen','PropSeen.Samples')])
total <- unique(total)
write.table(total, file = 'results/splice-events-summary/splice-events-rnaseq.txt', quote = F, sep = "\t", row.names = F)
