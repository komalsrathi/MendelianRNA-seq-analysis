# Author: Komal S. Rathi
# Date: 03/08/2019
# Function: script to further investigate 27571P

library(data.table)
setwd('~/Projects/DGD_Mendelian_RNASeq/')
dat <- fread('data/variant_filtering/rawdata/vardict/Sample_1__27571P.vardict-annotated-rnaedit-annotated-gemini-hgmdannotated.maf')
dat <- dat[which(dat$Hugo_Symbol == "NIPBL"),]
# vardict <- dat[grep('Splice|Intron', dat$Variant_Classification, ignore.case = T),c('Start_Position','End_Position','HGVSc','variant_qual','Variant_Classification','EXON','INTRON','RefSeq')]
vardict <- dat[,c('Start_Position','End_Position','HGVSc','variant_qual','Variant_Classification','EXON','INTRON','RefSeq')]
vardict$caller <- 'Vardict'

dat <- fread('data/variant_filtering/rawdata/gatk4/27571P-gatk-haplotype-annotated-hgmdannotated.maf')
dat <- dat[which(dat$Hugo_Symbol == "NIPBL"),]
gatk <- dat[,c('Start_Position','End_Position','HGVSc','variant_qual','Variant_Classification','EXON','INTRON','RefSeq')]
gatk$caller <- 'GATK'

dat <- fread('data/variant_filtering/rawdata/strelka/27571P.variants-hgmdannotated.maf')
dat <- dat[which(dat$Hugo_Symbol == "NIPBL"),]
strelka <- dat[,c('Start_Position','End_Position','HGVSc','variant_qual','Variant_Classification','EXON','INTRON','RefSeq')]
strelka$caller <- 'Strelka'

# these are all splice/intronic variants I could find in NIPBL by all three variant callers
total <- rbind(vardict, gatk, strelka)

total$splice_coords <- ifelse(total$Start_Position >= 36955739 & total$End_Position <= 36961586, 'Within_splice_coords','Outside')
total <- total[which(total$splice_coords == "Within_splice_coords"),]
write.table(total, file = 'data/27571P-splice-and-intronic-variants-in-NIPBL.txt', quote = F, sep = "\t", row.names = F)

# look up these variants in other samples
library(data.table)

# gatk
lf <- list.files(path = 'data/variant_filtering/rawdata/gatk4/', pattern = '*.maf', full.names = T)
x <- read.delim('data/27571P-splice-and-intronic-variants-in-NIPBL.txt', stringsAsFactors = F)
x <- unique(x$HGVSc)
for(i in 1:length(lf)){
  print(lf[i])
  dat <- fread(lf[i])
  dat <- dat[which(dat$Hugo_Symbol == "NIPBL"),]
  print(dat[which(dat$HGVSc %in% x),'HGVSc'])
}

# c.358+1477G>A in CDL-217-05P

# vardict
lf <- list.files(path = 'data/variant_filtering/rawdata/vardict/', pattern = '*.maf', full.names = T)
x <- read.delim('data/27571P-splice-and-intronic-variants-in-NIPBL.txt', stringsAsFactors = F)
x <- unique(x$HGVSc)
for(i in 1:length(lf)){
  print(lf[i])
  dat <- fread(lf[i])
  dat <- dat[which(dat$Hugo_Symbol == "NIPBL"),]
  print(dat[which(dat$HGVSc %in% x),'HGVSc'])
}

# c.358+1477G>A in CDL-217-05P

# strelka
lf <- list.files(path = 'data/variant_filtering/rawdata/strelka/', pattern = '*.maf', full.names = T)
x <- read.delim('data/27571P-splice-and-intronic-variants-in-NIPBL.txt', stringsAsFactors = F)
x <- unique(x$HGVSc)
for(i in 1:length(lf)){
  print(lf[i])
  dat <- fread(lf[i])
  dat <- dat[which(dat$Hugo_Symbol == "NIPBL"),]
  print(dat[which(dat$HGVSc %in% x),'HGVSc'])
}

# CDL-069-99P and CDL-123-01P has all three c.231-* variants
# c.327G>A is also present in CDL-515-09P
# c.358+1477G>A is also present in CDL-219-P