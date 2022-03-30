#Purpose:
#Identify if proteins that are HI and OE are enriched for stress granule formation
#In Yeast
#Here, the background we are using is the proteome

library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)


source('scripts/binning_functions.R')
source('scripts/enrichment_summary_functions.R')

master <- read.delim('data/merged_tables/yeast/actual_master.tsv',
                     sep = "\t")
nsg <- read.delim(file = 'data/phase_sep/Erich_nsg/nsg_yeast.txt', header = FALSE)
# log2+1 abundance
master$abundance <- log2(master$abundance+1)

### 1) Get Stress Granule Proteins

stress_granule_proteins <- filter(master, uniprot_id %in% nsg$V1) %>%
  mutate(stress_granule = as.integer(1))
not_stress_granule_proteins <- filter(master, !uniprot_id %in% nsg$V1) %>%
  mutate(stress_granule = as.integer(0))

### Append to master if the protein is a stress granule or not
# if it is, then the 'stress_granule' column will be == 1

master <- rbind(stress_granule_proteins, not_stress_granule_proteins)



#~~~~~~~~~ Global Parameters ~~~~~~~~~~~~

HI_col <- "HI"
OE_col <- "OE"

enrichment_column <-  'stress_granule'
alpha = 0.01


# #remove NAs
# master <- filter(master,
#                  !is.na(!!sym(abundance_col)) & !is.na(!!sym(disorder_col)))

### Get genes that are HI and OE

HI_OE <- filter(master, !!sym(HI_col) == 1 & !!sym(OE_col) == 1)

###### THeres only 20 genes that are HI and OE. Prob not high enough to do stats on?

