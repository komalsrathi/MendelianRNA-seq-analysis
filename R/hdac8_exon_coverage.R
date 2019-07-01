# Author: Komal S. Rathi
# Date: 08/30/2018
# Function: HDAC8 exon coverage analysis (check for RNA degradation in CDL-418)

setwd('~/Projects/DGD_Mendelian_RNASeq')
library(ggplot2)
library(tidyr)
library(dplyr)

lf <- list.files(pattern = '*_coverage.txt', path = 'results/HDAC8_analysis/hdac8-exon-coverage/coverage_results', full.names = T)
for(i in 1:length(lf)){
  print(i)
  fname <- gsub('.*/|_coverage.txt','',lf[i])
  if(i == 1){
    dat <- read.delim(lf[i], stringsAsFactors = F, header = F)
    dat$Sample <- fname
  } else {
    x <- read.delim(lf[i], stringsAsFactors = F, header = F)
    x$Sample <- fname
    dat <- rbind(dat, x)
  }
}

dat$exon_length <- dat$V3-dat$V2
dat <- dat[,c(4,7,8,9)]
colnames(dat) <- c("id","read_count","Sample","exon_length")
dat$transcript <- gsub('_.*','',dat$id)
dat$exon_n <- gsub('.*_','',gsub('_ENSE.*','',dat$id))
dat$exon_name <- gsub('.*_','',dat$id)

# add lib sizes
lib.size <- read.delim('results/HDAC8_analysis/hdac8-exon-coverage/flagstat_results/flagstat_output.txt', header = F)
dat <- merge(dat, lib.size, by.x = 'Sample', by.y = 'V1')
colnames(dat)[8] <- "library_size"

dat$per.million.scaling.factor <- dat$library_size/10^6
dat$FPKM <- ((dat$read_count+1)/dat$per.million.scaling.factor)/dat$exon_length
res <- unique(dat[,c(1,5,6,7,8,3,4,9,10)])

# plot transcript coverage boxplot
ggplot(res, aes(x = transcript, 
                y = read_count,
                group = Sample)) +
  geom_point(data = res, aes(color = transcript)) + 
  geom_line() + 
  facet_wrap(~Sample) +
  theme_bw() + 
  theme(axis.text.x = element_blank()) +
  ylab('Raw Read Counts') + xlab('') +
  guides(color=guide_legend(title="Transcript_ID"))
ggsave(filename = 'results/HDAC8_analysis/hdac8-exon-coverage/HDAC8_transcripts.pdf', width = 18, height = 15)  

# exon coverage for individual transcripts
plots <- function(x){
  pc <- c('ENST00000439122','ENST00000373573','ENST00000373589','ENST00000373554','ENST00000373556','ENST00000373559')
  tns <- unique(x$transcript)
  if(tns %in% pc){
    filename <- paste0('results/HDAC8_analysis/hdac8-exon-coverage/',tns,'_pc.pdf')
  } else {
    filename <- paste0('results/HDAC8_analysis/hdac8-exon-coverage/',tns,'.pdf')
  }
  ggplot(x, aes(x = Sample, y = read_count)) +
    geom_bar(stat = 'identity') +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap(~exon_name+exon_n) +
    ggtitle(tns)
  ggsave(filename = filename, width = 18, height = 15)  
}
plyr::d_ply(res, 'transcript', .fun = plots)

# canonical
lf <- list.files(pattern = '*_coverage.txt', path = 'results/HDAC8_analysis/hdac8-exon-coverage/coverage_results_v3', full.names = T)
for(i in 1:length(lf)){
  print(i)
  fname <- gsub('.*/|_coverage.txt','',lf[i])
  if(i == 1){
    dat <- read.delim(lf[i], stringsAsFactors = F, header = F)
    dat$Sample <- fname
  } else {
    x <- read.delim(lf[i], stringsAsFactors = F, header = F)
    x$Sample <- fname
    dat <- rbind(dat, x)
  }
}

dat$exon_length <- dat$V3-dat$V2
dat <- dat[,c(4,7,8,9)]
colnames(dat) <- c("id","read_count","Sample","exon_length")
dat$transcript <- gsub('_.*','',dat$id)
dat$exon_n <- gsub('.*_','',gsub('_ENSE.*','',dat$id))
dat$exon_name <- gsub('.*_','',dat$id)

# add lib sizes
lib.size <- read.delim('results/HDAC8_analysis/hdac8-exon-coverage/flagstat_results/flagstat_output.txt', header = F)
lib.size$V1 <- gsub("Sample_1__","",lib.size$V1)
dat <- merge(dat, lib.size, by.x = 'Sample', by.y = 'V1')
colnames(dat)[8] <- "library_size"

dat$per.million.scaling.factor <- dat$library_size/10^6
dat$FPKM <- ((dat$read_count+1)/dat$per.million.scaling.factor)/dat$exon_length
res <- unique(dat[,c(1,5,6,7,8,3,4,9,10)])

# plot transcript coverage boxplot
ggplot(res, aes(x = exon_n, 
                y = read_count,
                group = Sample)) +
  geom_point(data = res, aes(color = transcript)) + 
  geom_line() + 
  facet_wrap(~Sample) +
  theme_bw() + 
  theme(axis.text.x = element_blank()) +
  ylab('Raw Read Counts') + xlab('') +
  guides(color=guide_legend(title="Transcript_ID"))
ggsave(filename = 'results/HDAC8_analysis/hdac8-exon-coverage/HDAC8_transcripts_canonical.pdf', width = 18, height = 15)  
