# Author: Komal S Rathi
# Date: 04/30/2019
# Function: split variant caller plots 
# Step 2 (after using Mahdi's pipeline filters)

setwd('~/Projects/DGD_Mendelian_RNASeq/')
source('R/pubTheme.R')
library(reshape2)
library(ggplot2)

# samples
negative <- 'CDL-219-P'
samples <- read.delim('data/expression_data/patient_metadata.txt', stringsAsFactors = F)
samples <- samples[which(samples$id != negative),"id"]

# final list of variants from excel sheet
final <- read.delim('data/variant_filtering/final_splicevariants_exons.txt', stringsAsFactors = F)
final <- unique(final[,c("Sample","HGVSc","Exons","Introns","Label","Type","Hugo_Symbol")])

# caller <- 'vardict'
caller <- 'gatk3_source'

filtering.plot <- function(caller){
  
  # split final file into two
  splice <- final[which(final$Type == "Splice"),]
  nonsplice <- final[which(final$Type != "Splice"),]
  
  # read maf file with all variants
  maf.big <- paste0('data/variant_filtering/',caller,'_filter_breakdown_variants.maf')
  maf.big <- read.delim(maf.big, stringsAsFactors = F)
  
  # read maf file for final filtered variants
  maf <- paste0('data/variant_filtering/',caller,'_filtered_variants.maf')
  maf <- read.delim(maf, stringsAsFactors = F)
  maf$EXON <- gsub('[/].*', '',maf$EXON)
  maf$INTRON <- gsub('[/].*', '',maf$INTRON)
  maf <- maf[which(maf$Tumor_Sample_Barcode != negative),] # remove negative sample
  
  # subset calls to caller
  if(caller == "gatk4"){
    name <- "GATK4"
  } else if(caller == "vardict"){
    name <- "Vardict"
  } else if(caller == "gatk3") {
    name <- "GATK3"
  } else if(caller == "gatk3_source") {
    name <- "GATK3-source"
  } else if(caller == "strelka") {
    name <- "Strelka"
  }
  title <- paste0('Identification of Mutations via ', name)
  
  # output files
  fname <- paste0('paper/variant_plots/',name,'-filtering-plot.pdf')
  new.maf <- paste0('paper/variant_plots/',name,'-filter-breakdown.maf')
  maf.varclass <- paste0('paper/variant_plots/',name,'-variant-classes.txt')
  rawdata <- paste0('paper/variant_plots/',name,'-rawdata.txt')
  
  # df should be all sample + hugo symbol + exon matches
  for(i in 1:length(samples)){
    print(i)
    print(samples[i])
    splice.tmp <- splice[which(splice$Sample %in% samples[i]),]
    exons <- splice.tmp$Exons
    introns <- splice.tmp$Introns
    hgnc <- splice.tmp$Hugo_Symbol
    print(paste0("Exons to search: ", exons))
    print(paste0("Introns to search:", introns))
    print(paste0("Hugo symbol to search:", hgnc))
    exons <- unlist(strsplit(exons, ','))
    introns <- unlist(strsplit(introns, ','))
    if(i == 1){
      total <- maf[which(maf$Tumor_Sample_Barcode %in% samples[i]),]
      print(total[which(total$EXON %in% exons | total$INTRON %in% introns),c('Hugo_Symbol','EXON','INTRON')])
      total <- total[which(total$Hugo_Symbol %in% hgnc),] # newly added
      total <- total[which(total$EXON %in% exons | total$INTRON %in% introns),]
      print(total[,c('Hugo_Symbol','EXON','INTRON')])
      # print(paste0("Exons left: ",unique(total$EXON)))
      # print(paste0("Introns left: ",unique(total$INTRON)))
    } else {
      dat <- maf[which(maf$Tumor_Sample_Barcode %in% samples[i]),]
      print("Wrong:")
      print(dat[which(dat$EXON %in% exons | dat$INTRON %in% introns),c('Hugo_Symbol','EXON','INTRON')])
      dat <- dat[which(dat$Hugo_Symbol %in% hgnc),] # newly added
      dat <- dat[which(dat$EXON %in% exons | dat$INTRON %in% introns),]
      print("Correct:")
      print(dat[,c('Hugo_Symbol','EXON','INTRON')])
      # print(paste0("Exons left: ", unique(dat$EXON)))
      # print(paste0("Introns left: ",unique(total$INTRON)))
      total <- rbind(total, dat)
    }
  }
  rm(dat)
  df <- total
  maf.big[which(maf.big$var_id %in% df$var_id),"F4"] <- "Y"
  write.table(maf.big, file = new.maf, quote = F, sep = "\t", row.names = F)
  df <- unique(df[,c('Tumor_Sample_Barcode','HGVSc')])
  
  # get overlaps between filters and final splice variants
  n.final <- df[which(df$Tumor_Sample_Barcode %in% splice$Sample & df$HGVSc %in% splice$HGVSc),]
  n.final <- n.final$Tumor_Sample_Barcode
  print(n.final)
  df <- plyr::count(df$Tumor_Sample_Barcode)
  remaining.samples <- setdiff(samples, df$x)
  df <- rbind(df, data.frame(x = remaining.samples, freq = 0))
  df$variable <- 'F4'
  df$value <- df$freq
  df$sample <- df$x 
  df$label <- df$value
  df$label <- ifelse(df$sample %in% n.final, paste0(df$label,"*"), df$label)
  # df$label <- ifelse(df$sample == "CDL-223-05P", paste0(df$label,"*"), df$label) # manually add CDL-223
  
  # filtered variants numbers
  vars.filtered <- paste0('data/variant_filtering/',caller,'_filtered_variants.txt')
  vars.filtered <- read.delim(vars.filtered, stringsAsFactors = F)
  vars.filtered$exonic <- NULL
  vars.filtered$intronic <- NULL
  vars.filtered <- melt(vars.filtered)
  vars.filtered <- vars.filtered[which(vars.filtered$sample != negative),]
  vars.filtered$label <- vars.filtered$value
  
  # for non-splice add * in F3
  non.splice.samples <- maf[which(maf$Tumor_Sample_Barcode %in% nonsplice$Sample & maf$HGVSc %in% nonsplice$HGVSc),'Tumor_Sample_Barcode']
  vars.filtered$label <- ifelse(vars.filtered$sample %in% non.splice.samples & vars.filtered$variable == "F3", paste0(vars.filtered$label,"*"), vars.filtered$label)
  
  # bind both
  df <- df[,colnames(vars.filtered)]
  vars.filtered <- rbind(vars.filtered, df)
  
  vars.filtered$annotation <- ifelse(vars.filtered$sample %in% c('95-0614-P',
                                                                 'CDL-022-P',
                                                                 'CDL-086-99P',
                                                                 'CDL-069-99P',
                                                                 'CDL-679-14P'), 'N', 'P')
  
  vars.filtered$sample <- paste0(vars.filtered$sample, ' (', vars.filtered$annotation, ')')
  vars.filtered$variable <- as.character(vars.filtered$variable)
  
  if(caller == "strelka"){
    scale.max <- max(vars.filtered$value) + 100
  } else {
    scale.max <- max(vars.filtered$value) + 10
  }
  write.table(vars.filtered, file = rawdata, quote = F, sep = "\t", row.names = F)
  
  # before plotting, reorder labels
  vars.filtered$sample <- sub('-[0-9]{0,2}[A-Z]{0,1} [(]',' (', vars.filtered$sample)
  vars.filtered <- vars.filtered[order(vars.filtered$annotation, decreasing = T),]
  vars.filtered$sample <- ifelse(vars.filtered$sample %in% c("CDL-223 (P)","27571P (P)"), paste0(vars.filtered$sample,"**"), vars.filtered$sample)
  vars.filtered$sample <- factor(vars.filtered$sample, levels = unique(vars.filtered$sample))
  
  p <- ggplot(vars.filtered, aes(x = variable, y = value, group = sample)) + 
    geom_line(linetype = "dashed") + 
    geom_point(aes(color = variable)) + 
    geom_text(aes(label = label, vjust = -0.3, hjust = 0), 
              position = position_dodge(width = 1), size = 3) +
    facet_wrap(~sample, nrow = 3) + theme_Publication3() +
    scale_color_discrete(name = "Filters",
                         breaks = c("F1","F2","F3","F4"),
                         labels = c("F1 = Variants within exons +/- 10bp", 
                                    "F2 = Variant Quality >= 30", 
                                    "F3 = gnomAD AF filter (all subpopulations)\n<0.5% HGMD DM/DM? variants\n<0.1% non-HGMD",
                                    "F4 = Flanking exons (+/- 1 exon) of known splice variants")) +
    xlab('Filtering Steps') + ylab('Number of Mutations') +
    ggtitle(title) + 
    scale_y_continuous(limits = c(0, scale.max)) 
  p
  ggsave(filename = fname, plot = p, width = 12, height = 8)
  
  # write down variant classification for F6
  maf <- plyr::count(maf, c("Tumor_Sample_Barcode","Variant_Classification"))
  maf <- dcast(maf, Tumor_Sample_Barcode~Variant_Classification, value.var = 'freq')
  maf[is.na(maf)] <- 0
  maf$Total <- rowSums(maf[,2:ncol(maf)])
  write.table(x = maf, file = maf.varclass, quote = F, sep = "\t", row.names = F)
}

filtering.plot(caller = 'gatk4')
filtering.plot(caller = 'gatk3') # gatk3.8 (bcbio)
filtering.plot(caller = 'gatk3_source') # gatk3.8 (source)
filtering.plot(caller = 'vardict')
filtering.plot(caller = 'strelka')

# this is to double check if all the variants were captured:
final <- read.delim('data/variant_filtering/final_splicevariants_exons.txt', stringsAsFactors = F)
hgvsc <- unique(final$HGVSc)
one <- read.delim('data/variant_filtering/gatk_filtered_variants.maf', stringsAsFactors = F)
two <- read.delim('data/variant_filtering/vardict_filtered_variants.maf', stringsAsFactors = F)
three <- read.delim('data/variant_filtering/strelka_filtered_variants.maf', stringsAsFactors = F)
dat <- rbind(one, two, three)
dat[which(dat$HGVSc %in% hgvsc),'HGVSc']
dat <- dat[which(dat$HGVSc %in% hgvsc),]
setdiff(hgvsc, dat$HGVSc)
