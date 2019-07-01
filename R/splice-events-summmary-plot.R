# Author: Komal S. Rathi
# Date: 04/30/2019
# Function: Splice events summary (dot-line) plot

setwd('~/Projects/DGD_Mendelian_RNASeq/')
source('R/pubTheme.R')
library(reshape2)
library(dplyr)
library(tidyr)
library(ggplot2)

# total number of splice variants called
negative <- 'CDL-219-P'
patients.ct <- read.delim('results/splice-events-summary/patient_splice_events.txt')
patients.ct <- patients.ct[,c(1,3)]
patients.ct <- patients.ct[which(patients.ct$Samples != negative),]
colnames(patients.ct) <- c("Sample","freq")
patients.ct$filter <- 'F1'

# step 1
filt.reads <- read.delim('results/splice-events-patient/FilteredJunctions_ReadSupport10_FilterReads5.txt')
filt.reads <- filt.reads[-which(filt.reads$Sample %in% c('-',negative)),]
filt.reads <- filt.reads[which(filt.reads$NTimesSeen >= 10), ]
step1 <- plyr::count(filt.reads$Sample)
step1$filter <- 'F2'
colnames(step1)[1] <- 'Sample'

# step 2
only1sample <- read.delim('results/splice-events-patient/FilteredJunctions_ReadSupport10_Only1Sample.txt')
only1sample <- unique(only1sample[,c('Sample','Pos')])
step2 <- plyr::count(only1sample$Sample)
to.add <- setdiff(patients.ct$Sample, step2$x)
to.add <- data.frame(x = to.add, freq = 0)
step2 <- rbind(step2, to.add)
colnames(step2)[1] <- 'Sample'
step2$filter <- 'F3'

# step 3
# step3 <- filt.reads[which(filt.reads$NSamplesSeen == 1),]
# step3 <- step3[which(step3$NTimesSeen >= 10),]
# step3 <- step3[-which(step3$Pos %in% only1sample$Pos),]
# step3 <- unique(step3[,c("Sample","Pos")])
# step3 <- plyr::count(step3$Sample)
# to.add <- setdiff(patients.ct$Sample, step3$x)
# to.add <- data.frame(x = to.add, freq = 0)
# step3 <- rbind(step3, to.add)
# colnames(step3)[1] <- 'Sample'
# step3$filter <- 'F4'

# step 4: number final
final <- read.delim('results/splice-events-final-summary.txt')
final$freq <- ifelse(final$Decision == "Positive", 1, 0)
final$Decision <- NULL
final$filter <- 'F4'

# step 5: dna
dna <- final
dna$filter <- 'F5'

# plot
total <- rbind(patients.ct, step1, step2, final, dna)
lev <- unique(total$filter)
total$filter <- factor(total$filter, levels = lev)

# add annotation to samples
total$annotation <- ifelse(total$Sample %in% c('95-0614-P',
                                               'CDL-022-P',
                                               'CDL-086-99P',
                                               'CDL-069-99P',
                                               'CDL-679-14P'), 'N', 'P')
total$Sample <- paste0(total$Sample, ' (', total$annotation, ')')

# order by sample type
total$Sample <- sub('-[0-9]{0,2}[A-Z]{0,1} [(]',' (', total$Sample)
total <- total[order(total$annotation, decreasing = T),]
total$Sample <- factor(total$Sample, levels = unique(total$Sample))

# plot for splice events
p <- ggplot(total, aes(x = filter, y = freq, group = Sample)) + 
  geom_point(aes(color = filter)) + geom_line(linetype = "dotted") + 
  geom_text(aes(label = freq, vjust = -0.3, hjust = 0), 
            position = position_dodge(width = 1)) +
  facet_wrap(~Sample, nrow = 3) + theme_Publication3() +
  scale_color_discrete(name="Filters",
                       breaks = c("F1","F2","F3","F4","F5"),
                       labels=c("F1 = Total splice events", 
                                'F2 = Total read support >= 10', 
                                'F3 = Observed in only 1 Patient',
                                'F4 = Non-overlapping in GTEx Whole Blood and EBV-transformed lymphocytes samples',
                                'F5 = Most likely true genomic variant')) +
  xlab('Filtering Steps') + ylab('Number of Splice Events') +
  ggtitle('Identification of deleterious splice events') + 
  scale_y_continuous(limits = c(0, 1000)) 
p
ggsave(filename = 'paper/splice-events-filtering-plot.pdf', plot = p, width = 12, height = 8)

# events <- data.frame(filter = c("F3","F5"), n = c(3,5))
# geom_vline(data = events, aes(xintercept = n), linetype="dotted", color = 'red')
