# Author: Komal S. Rathi
# Date: 06/21/2018
# Function: Rename strelka output (on server)
# replace "variants.vcf.gz" to "samplename.vcf.gz"

setwd('/mnt/isilon/cbmi/variome/rathik/mendelian_rnaseq/strelka_out')
dirs <- list.dirs(path = '/mnt/isilon/cbmi/variome/rathik/mendelian_rnaseq/strelka_out', recursive = F, full.names = T)

# replace all variants.vcf.gz
for(i in 1:length(dirs)){
  print(i)
  full.path <- paste0(dirs[i],'/results/variants/variants.vcf.gz')
  sample.name <- gsub('.*/','',dirs[i])
  replace.path <- paste0(dirs[i],'/results/variants/',sample.name,'.variants.vcf.gz')
  print(paste0("Original: ", full.path))
  print(paste0("Replacement: ", replace.path))
  file.rename(from = full.path, to = replace.path)
}

for(i in 1:length(dirs)){
  print(i)
  full.path <- paste0(dirs[i],'/results/variants/genome.S1.vcf.gz')
  sample.name <- gsub('.*/','',dirs[i])
  replace.path <- paste0(dirs[i],'/results/variants/',sample.name,'.genome.S1.vcf.gz')
  print(paste0("Original: ", full.path))
  print(paste0("Replacement: ", replace.path))
  file.rename(from = full.path, to = replace.path)
}

total <- 'x'
for(i in 1:length(dirs)){
  print(i)
  sample.name <- gsub('.*/','',dirs[i])
  path1 <- paste0(dirs[i],'/results/variants/',sample.name,'.variants.vcf.gz')
  if(file.exists(path1)){
    total <- c(total, path1)
  }
}
write.table(total, file = '/mnt/isilon/cbmi/variome/rathik/mendelian_rnaseq/strelka_out/filepaths.txt', quote = F, sep = "\t", row.names = F)