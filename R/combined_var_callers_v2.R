# Author: Komal S Rathi
# Date: 04/30/2019
# Function: Combined Variant caller plots 
# Step 3 (after individual plots)
# (GATK3.8, GATK4 and Vardict)

setwd('~/Projects/DGD_Mendelian_RNASeq/')
source('R/pubTheme.R')
library(ggplot2)

# combined barplot
strelka <- read.delim('paper/variant_plots/Strelka-rawdata.txt', stringsAsFactors = F)
strelka$caller <- 'Strelka'
vardict <- read.delim('paper/variant_plots/Vardict-rawdata.txt', stringsAsFactors = F)
vardict$caller <- 'Vardict'
gatk <- read.delim('paper/variant_plots/GATK4-rawdata.txt', stringsAsFactors = F)
gatk$caller <- 'GATK4'
gatk3 <- read.delim('paper/variant_plots/GATK3-source-rawdata.txt', stringsAsFactors = F)
gatk3$caller <- "GATK3.8"
total <- rbind(gatk3, gatk, vardict, strelka)

# order by sample type
total$sample <- sub('-[0-9]{0,2}[A-Z]{0,1} [(]',' (', total$sample)
total <- total[order(total$annotation, decreasing = T),]
total$sample <- ifelse(total$sample %in% c("CDL-223 (P)","27571P (P)"), paste0(total$sample,"**"), total$sample)
total$sample <- factor(total$sample, levels = unique(total$sample))

# order by variant callers
total$caller <- factor(total$caller, levels = c('Strelka','GATK3.8','GATK4','Vardict'))
total <- total[which(total$caller != "Strelka"),]

# plot 
# scale.max <- max(log2(total$value + 1)) + 2
scale.max <- max(total$value) + 5
p <- ggplot(total, aes(x = variable, y = value, group = caller)) + 
  geom_bar(position = position_dodge(width=0.8), stat="identity", 
           width=0.8, aes(fill = variable), color = 'black', size = 0.2) +
  geom_text(aes(label = label, color= caller, vjust = -0.5, hjust = 0.5),
            position = position_dodge(width = 0.8), size = 2) +
  facet_wrap(~sample, nrow = 3) + theme_Publication4() + 
  scale_color_manual(name = "Caller",
                     values = c("GATK3.8" = "black", 
                                "Vardict" = "blue3", 
                                "GATK4" = "red")) +
  scale_fill_discrete(name = "Filters",
                      breaks = c("F1","F2","F3","F4"),
                      labels = c("F1 = Variants within exons +/- 10bp", 
                                 "F2 = Variant Quality >= 30", 
                                 "F3 = gnomAD AF filter (all subpopulations)\n       <0.5% HGMD DM/DM? variants\n       <0.1% non-HGMD",
                                 "F4 = Flanking exons (+/- 1 exon) of known splice variants")) +
  xlab('Filtering Steps') + ylab('Number of Variants (log2)') +
  ggtitle('Variant Calling Filtering Pipeline') + 
  scale_y_continuous(limits = c(0, scale.max)) 
p
ggsave(filename = 'paper/variant_plots/combined_variant_plot_v2.pdf', plot = p, width = 10, height = 6)

