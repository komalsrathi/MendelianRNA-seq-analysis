####################################################################
# Author: Komal S Rathi
# Date: 01/31/2019
# Function: script to filter variants from GATK, Vardict and Strelka
# Mahdi's pipeline filters
# 1. ROI filter: exons +/- 10
# 2. Qual by depth: 5 for low GQ/low DP variants
# 3. Population Filters:
#   HGMD variants 1%, synonymous variants 0.1%, other variants 0.5%
# how many remain - save to another file
# Step 1 (after annotating maf with HGMD)
####################################################################

library(data.table)
library(hutils)
library(tidyr)
library(GenomicRanges)
library(reshape2)
library(dplyr)

# setwd('/mnt/isilon/cbmi/variome/rathik/mendelian_rnaseq')
setwd('~/Projects/DGD_Mendelian_RNASeq/')

# use exon file from biomart output
exon.dat <- read.delim('data/gencode.v19.cdl_canonical_transcripts.v7.patched_contigs.exons.txt', stringsAsFactors = F)
exon.dat <- exon.dat %>% group_by(gene_symbol) %>%
  mutate(exon_start = exon_start - 10,
         exon_end = exon_end + 10) %>%
  unique() %>% dplyr::select(gene_symbol, chrom, exon_start, exon_end) %>%
  as.data.frame()

# final exons/intron list
# introns <- read.delim('data/variant_filtering/final_splicevariants_exons.txt', stringsAsFactors = F)
# introns <- unique(introns[,c("Sample","HGVSc","Introns","Label","Hugo_Symbol")])

# vars to test (for testing)
# vars.to.test <- c("c.64","c.359","c.4561")
# vars.to.test <- paste(vars.to.test, collapse = "|")

# folder
folder <- 'data/variant_filtering/rawdata/gatk3-source/'

filter.out <- function(folder){
  lf <- list.files(path = folder, pattern = '*.maf', full.names = TRUE)
  
  for(i in 1:length(lf)){
    print(paste0("Sample no.: ",i))
    n <- gsub('.*/|-hgmdannotated.maf|-gatk-haplotype-annotated-hgmdannotated.maf|.variants-hgmdannotated.maf|.vardict-annotated-rnaedit-annotated-gemini-hgmdannotated.maf|Sample_1__','',lf[i])
    print(paste0("Sample: ",n))
    dat <- data.table::fread(lf[i], verbose = FALSE)
    dat <- as.data.frame(dat)
    dat <- dat[which(dat$variant_qual != "."),]
    dat$Tumor_Sample_Barcode <- n # add Tumor Sample Barcode
    dat$var_id <- paste0(dat$Tumor_Sample_Barcode,'_', rownames(dat))
    
    # ROI filter:
    # identify all exonic variants
    # add sequential identifiers
    exonic.vars <- dat
    exonic.vars$id <- seq(1:nrow(exonic.vars))
    
    # only keep positions that are within the exons in the exon file
    subject <- with(exon.dat, GRanges(chrom, IRanges(start = exon_start, end = exon_end, names = gene_symbol)))
    query <- with(exonic.vars, GRanges(Chromosome, IRanges(start = Start_Position, end = End_Position, names = id)))
    
    # find overlaps and subset maf 
    res <- findOverlaps(query = query, subject = subject, type = "within")
    res.df <- data.frame(exonic.vars[queryHits(res),], exon.dat[subjectHits(res),])
    exonic.vars <- exonic.vars[which(exonic.vars$id %in% res.df$id),]
    exonic.vars$id <- NULL
    print(paste0("Dimensions of exonic vars: ", nrow(exonic.vars)))
    dat <- exonic.vars
    maf <- dat # this is for writing out full maf
    maf[which(maf$var_id %in% dat$var_id),"F1"] <- "Y"
    s0 <- nrow(dat)
    print(s0)
    
    # Quality Filters:
    # quality by depth >= 2
    dat$variant_qual <- as.numeric(dat$variant_qual)
    dat <- dat[which(dat$variant_qual >= 30),]
    maf[which(maf$var_id %in% dat$var_id),"F2"] <- "Y"
    s1 <- nrow(dat) 
    print(s1)
    
    # Population Filters
    # if HGMD annotated, AF <= 0.01 else Syn variants <= 0.001 and others <= 0.005
    # replace gnomAD NAs with 0s
    dat[,grep('gnomAD_.*AF$', colnames(dat))][is.na(dat[,grep('gnomAD_.*AF$', colnames(dat))])] <- 0 
    dat$AF_filter <- FALSE
    dat$gnomAD_max_AF <- apply(dat[,grep('gnomAD_[A-Z]{3}_AF',colnames(dat))], 1, max)
    print(summary(dat$gnomAD_max_AF))
    dat$AF_filter <- ifelse(is.na(dat$CLASS),
                            ifelse(dat$gnomAD_max_AF <= 0.001, TRUE, FALSE), 
                            ifelse(dat$gnomAD_max_AF <= 0.005, TRUE, FALSE))
    dat <- dat[which(dat$AF_filter == TRUE),]
    maf[which(maf$var_id %in% dat$var_id),"F3"] <- "Y"
    s2 <- nrow(dat)
    print(s2)
    
    # add sample name to tumor_sample_barcode
    t <- data.frame(sample = n, exonic = nrow(exonic.vars), F1 = s0, F2 = s1, F3 = s2)
    
    if(i == 1){
      total1 <- t
      total2 <- dat
      total3 <- maf
    } else {
      total1 <- rbind(total1, t)
      total2 <- rbind(total2, dat)
      total3 <- rbind(total3, maf)
    }
  }
  
  # return results
  return(list(total1, total2, total3))
  
}

# exonic pipeline
vardict.total <- filter.out(folder = 'data/variant_filtering/rawdata/vardict/')
write.table(vardict.total[[1]], file = 'data/variant_filtering/vardict_filtered_variants.txt', quote = F, sep = "\t", row.names = F)
write.table(vardict.total[[2]], file = 'data/variant_filtering/vardict_filtered_variants.maf', quote = F, sep = "\t", row.names = F)
write.table(vardict.total[[3]], file = 'data/variant_filtering/vardict_filter_breakdown_variants.maf', quote = F, sep = "\t", row.names = F)

gatk4.total <- filter.out(folder = 'data/variant_filtering/rawdata/gatk4/')
write.table(gatk4.total[[1]], file = 'data/variant_filtering/gatk4_filtered_variants.txt', quote = F, sep = "\t", row.names = F)
write.table(gatk4.total[[2]], file = 'data/variant_filtering/gatk4_filtered_variants.maf', quote = F, sep = "\t", row.names = F)
write.table(gatk4.total[[3]], file = 'data/variant_filtering/gatk4_filter_breakdown_variants.maf', quote = F, sep = "\t", row.names = F)

strelka.total <- filter.out(folder = 'data/variant_filtering/rawdata/strelka/')
write.table(strelka.total[[1]], file = 'data/variant_filtering/strelka_filtered_variants.txt', quote = F, sep = "\t", row.names = F)
write.table(strelka.total[[2]], file = 'data/variant_filtering/strelka_filtered_variants.maf', quote = F, sep = "\t", row.names = F)
write.table(strelka.total[[3]], file = 'data/variant_filtering/strelka_filter_breakdown_variants.maf', quote = F, sep = "\t", row.names = F)

# gatk 3.8 (bcbio)
gatk3.total <- filter.out(folder = 'data/variant_filtering/rawdata/gatk3/')
write.table(gatk3.total[[1]], file = 'data/variant_filtering/gatk3_filtered_variants.txt', quote = F, sep = "\t", row.names = F)
write.table(gatk3.total[[2]], file = 'data/variant_filtering/gatk3_filtered_variants.maf', quote = F, sep = "\t", row.names = F)
write.table(gatk3.total[[3]], file = 'data/variant_filtering/gatk3_filter_breakdown_variants.maf', quote = F, sep = "\t", row.names = F)

# gatk 3.8 (source)
gatk3.total <- filter.out(folder = 'data/variant_filtering/rawdata/gatk3-source/')
write.table(gatk3.total[[1]], file = 'data/variant_filtering/gatk3_source_filtered_variants.txt', quote = F, sep = "\t", row.names = F)
write.table(gatk3.total[[2]], file = 'data/variant_filtering/gatk3_source_filtered_variants.maf', quote = F, sep = "\t", row.names = F)
write.table(gatk3.total[[3]], file = 'data/variant_filtering/gatk3_source_filter_breakdown_variants.maf', quote = F, sep = "\t", row.names = F)
