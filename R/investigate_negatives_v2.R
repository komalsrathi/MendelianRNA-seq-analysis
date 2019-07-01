# Author: Komal S. Rathi
# Date: 09/10/2018
# Function: Investigate negative samples that have splice variants (v2)

library(reshape2)
library(ggplot2)
library(data.table)

setwd('~/Projects/DGD_Mendelian_RNASeq/results/mafs/strelka/')
# 95-0614 exon gain between Exon 21-22 chr5:37010041-37014938
# gatk - nothing reported
# dat <- fread('results/mafs/gatk/95-0614-P-gatk-haplotype-annotated.maf')
# res <- dat[which(dat$Hugo_Symbol == "NIPBL" & dat$Start_Position >= 37010041 & dat$Start_Position <= 37014938),]

# vardict - nothing reported
# dat <- fread('results/mafs/vardict/Sample_1__95-0614-P.vardict-annotated-rnaedit-annotated-gemini.maf')
# res <- dat[which(dat$Hugo_Symbol == "NIPBL" & dat$Start_Position >= 37010041 & dat$Start_Position <= 37014938),]

# strelka
dat <- fread('data/95-0614-P.variants.maf')
res <- dat[which(dat$Hugo_Symbol == "NIPBL" & dat$Start_Position >= 37010041 & dat$Start_Position <= 37014938),]
res <- res[,c("Hugo_Symbol","Chromosome","Start_Position","Variant_Classification", "Variant_Type", "dbSNP_RS", "HGVSc", "HGVSp_Short", "Exon_Number", "t_depth", "t_ref_count","t_alt_count","FILTER")]
print(nrow(res))
# write.table(res, file = 'results/95-0614-P_variants.txt', quote = F, sep = "\t", row.names = F)

# test which variants are present in other samples
lf <- list.files('data/', pattern = '*.maf', full.names = T)
for(i in 1:length(lf)){
  print(lf[i])
  strelka.2 <- fread(lf[i], verbose = FALSE)
  strelka.2 <- strelka.2[which(strelka.2$Hugo_Symbol == "NIPBL" & strelka.2$Start_Position >= 37010041 & strelka.2$Start_Position <= 37014938),]
  print(intersect(res$Start_Position, strelka.2$Start_Position))
}

# CDL-679 mutation Exon 33 and 34 skipping (5:37027514-37044449)
strelka <- fread('data/CDL-679-14P.variants.maf')
strelka.res <- strelka[which(strelka$Hugo_Symbol == "NIPBL" & strelka$Start_Position >= 37027514 & strelka$Start_Position <= 37044449),]
strelka.res <- strelka.res[,c("Hugo_Symbol","Chromosome","Start_Position","Variant_Classification", "Variant_Type", "dbSNP_RS", "HGVSc", "HGVSp_Short", "Exon_Number", "t_depth", "t_ref_count","t_alt_count","FILTER")]
print(nrow(strelka.res))
# write.table(strelka.res, file = 'results/CDL-679-14P_variants.txt', quote = F, sep = "\t", row.names = F)

# gatk <- fread('results/mafs/gatk/CDL-679-14P-gatk-haplotype-annotated.maf')
# gatk.res <- gatk[which(gatk$Hugo_Symbol == "NIPBL" & gatk$Start_Position >= 37027514 & gatk$Start_Position <= 37044449),]
# gatk.res <- gatk.res[,c("Hugo_Symbol","Chromosome","Start_Position","Variant_Classification", "Variant_Type", "dbSNP_RS", "HGVSc", "HGVSp_Short", "Exon_Number", "t_depth", "t_ref_count","t_alt_count","FILTER")]

# vardict <- fread('results/mafs/vardict/Sample_1__CDL-679-14P.vardict-annotated-rnaedit-annotated-gemini.maf')
# vardict.res <- vardict[which(vardict$Hugo_Symbol == "NIPBL" & vardict$Start_Position >= 37027514 & vardict$Start_Position <= 37044449),]
# vardict.res <- vardict.res[,c("Hugo_Symbol","Chromosome","Start_Position","Variant_Classification", "Variant_Type", "dbSNP_RS", "HGVSc", "HGVSp_Short", "Exon_Number", "t_depth", "t_ref_count","t_alt_count","FILTER")]

# test which variants are present in other samples
lf <- list.files('.', pattern = '*.maf', full.names = T)
for(i in 1:length(lf)){
  print(lf[i])
  strelka.2 <- fread(lf[i], verbose = FALSE)
  strelka.2 <- strelka.2[which(strelka.2$Hugo_Symbol == "NIPBL" & strelka.2$Start_Position >= 37027514 & strelka.2$Start_Position <= 37044449),]
  print(intersect(strelka.res$Start_Position, strelka.2$Start_Position))
}
