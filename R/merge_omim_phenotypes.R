# Author: Komal S. Rathi
# Date: 04/24/2019
# Function: 
# Merge OMIM and Phenotypic Series (from omim api)

# add omim info to the table (4141 genes)
omim <- data.table::fread('data/genemap2.txt') # from omim api
omim <- omim %>% mutate(Phenotypes = strsplit(as.character(Phenotypes),";")) %>%
  unnest(Phenotypes) %>%
  unique() %>%
  as.data.frame()
omim$Phenotypes <- gsub('[(].*','',omim$Phenotypes)
omim$Pheno.MIM.number <- omim$Phenotypes
omim$Pheno.MIM.number <- gsub('.*, ','',omim$Pheno.MIM.number)
omim$Pheno.MIM.number <- gsub(' $','',omim$Pheno.MIM.number)
omim$Pheno.MIM.number <- apply(omim, 1, FUN = function(x) ifelse(length(grep('[A-Za-z]',x[15]))>0,'',x[15]))
omim$Phenotypes <- gsub(', [0-9]{6}','',omim$Phenotypes)
omim$Phenotypes <- trimws(omim$Phenotypes)
omim <- unique(omim[,c('Ensembl Gene ID', 'Phenotypes','Pheno.MIM.number','Mim Number')])
omim <- omim[which(omim$`Ensembl Gene ID` != "" & omim$Phenotypes != ""),]

# add phenotypic series
pheno <- data.table::fread('data/phenotypicSeries_v2.txt', skip = 2) # omim api
pheno$set <- apply(pheno, 1, FUN = function(x) ifelse(length(grep('[A-Za-z]', x[2]))>0,'Superset','Subset'))
superset <- pheno[which(pheno$set == "Superset"),]
subset <- pheno[which(pheno$set == "Subset"),]
colnames(superset)[2] <- "Phenotypic_Superset"
pheno <- merge(superset[,1:2], subset[,1:2], by = 'Phenotypic Series Number')

# merge omim and phenotypic series
add.pheno <- merge(pheno, omim, by.x = 'MIM Number', by.y = 'Pheno.MIM.number', all.y = TRUE)
add.pheno$Phenotypic_Superset <- ifelse(is.na(add.pheno$Phenotypic_Superset), add.pheno$Phenotypes, add.pheno$Phenotypic_Superset)
add.pheno <- add.pheno %>% group_by(Phenotypic_Superset) %>%
  unique() %>%
  mutate(total.genes = n()) %>%
  unique() %>%
  as.data.frame()
write.table(add.pheno, file = 'data/genemap_phenotypicSeries_merged.txt', quote = F, sep = "\t", row.names = F)


