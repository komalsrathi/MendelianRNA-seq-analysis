# Author: Komal S. Rathi
# Date: 06/12/2018
# Function: Barplots for CdLS genes (only HDAC8, SMC1A and SMC3) RNA expression
# check for RNA degradation

library(RDiseaseXpress)
library(dplyr)
library(ggplot2)

setwd('~/Projects/DGD_Mendelian_RNASeq/')
source('R/pubTheme.R')

# gene expression data
genes <- c("SMC1A","HDAC8","SMC3")
dat <- getDataAnnotationByGeneSymbol(myGeneSymbols = genes, myStudy = "GTEx", myNorms = "rsem")
dat <- unique(dat[,c("gene_symbol","gene_id","data.sample_id","data.rsem.fpkm")])
annot <- unique(dat[,c('gene_id','gene_symbol')])

# meta
setwd('~/Projects/DGD_Mendelian_RNASeq/')
meta <- read.delim('data/gtex_data/GTEx_metadata_593_RNAseq.txt')
dat.sub <- merge(dat, meta, by.x = 'data.sample_id', by.y = 'Run')
dat.sub <- dat.sub %>% 
  group_by(gene_id, body_site) %>% 
  summarise(median = median(data.rsem.fpkm)) %>%
  as.data.frame()
dat.sub <- dat.sub[,c('gene_id','median','body_site')]
colnames(dat.sub) <- c("gene_id","FPKM",'Sample')

# barplot
p <- readRDS('results/expression/CDLS_patient_expr_hg38.RDS')
p <- p[which(p$gene_id %in% dat.sub$gene_id),c('gene_id','FPKM','Sample')]
total <- rbind(dat.sub, p)
total$Sample <- as.character(total$Sample)
total$Sample[total$Sample == "Cells - EBV-transformed lymphocytes"] <- "EBV-transformed"
total <- merge(total, annot, by.x = 'gene_id')

# expression plot
q <- ggplot(total, aes(Sample, FPKM)) + 
  geom_bar(stat = 'identity') + theme_Publication() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  facet_wrap(~gene_symbol, ncol = 1, scales = "free") +
  ggtitle("Gene Expression")
ggsave(q, filename = 'results/expression/gene_expr.pdf', device = 'pdf', width = 8, height = 10)  

# transcript level changes
p <- readRDS('results/expression/CDLS_patient_expr_transcripts_hg38.RDS')
p <- p[which(p$gene_id %in% annot$gene_id),]

for(i in 1:nrow(annot)){
  gi <- annot[i,1]
  gs <- annot[i,2]
  tmp <- p[which(p$gene_id %in% gi),]
  fname <- paste0('results/expression/',gs,'_transcripts.pdf')
  t <- ggplot(tmp, aes(transcript_id, TPM)) + 
    geom_bar(stat = 'identity') + theme_Publication() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
    facet_wrap(~Sample, ncol = 4, scales = "free") +
    ggtitle(paste0(gs," (TPM)")) + ylab("") +
    theme(plot.margin = unit(c(.5,.5,.5,2), "cm"))
  if(i == 3){
    ggsave(filename = fname, plot = t, device = 'pdf', width = 25, height = 10)
  } else {
    ggsave(filename = fname, plot = t, device = 'pdf', width = 15, height = 10)
  }
}
