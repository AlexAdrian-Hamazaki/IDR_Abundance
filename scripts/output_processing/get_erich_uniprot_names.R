library(tidyr)

#run this script to get binned_genes
source(file = 'scripts/yeast/analysis/HI/disopred3/bin_enrichment_of_HI.R')



get_uniprot_names <- function(ttt) {
  # for use in lapply loop
  # Gets uniprot ids for each bin.

  return(ttt[,1])
}


uniprot_names <- lapply(binned_genes, get_uniprot_names)

saveRDS(uniprot_names, file = 'data/for_erich/Yeast5by5')

a <- readRDS(file = 'data/for_erich/Human3by3')
b <- readRDS(file = 'data/for_erich/Human4by4')
c <- readRDS(file = 'data/for_erich/Human5by5')
d <- readRDS(file = 'data/for_erich/Yeast3by3')
e <- readRDS(file = 'data/for_erich/Yeast4by4')
f <- readRDS(file = 'data/for_erich/Yeast5by5')