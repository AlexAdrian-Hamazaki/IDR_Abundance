#Purpose:

#Combine human_total and HI bins and save it as human_master.R

library(dplyr)
library(tidyr)
library(stringr)

human_total <- read.delim('data/merged_tables/human/human_total_data.tsv')
HI_bins <- read.delim(file = 'data/merged_tables/human/human_HI_bins.tsv')
cancer <- read.delim(file = "data/dose/human/Census_allMon Oct 25 22 33 06 2021.tsv",
                        sep = "\t")

human_master <- left_join(human_total, HI_bins, by = "uniprot_id")

cancer <- cancer %>%
  select("Gene.Symbol", "Role.in.Cancer", "Entrez.GeneId") %>%
  rename("uniprot_symbol" = "Gene.Symbol")

sum(str_detect(cancer$Role.in.Cancer, 'oncogene'))

human_master <- left_join(human_master, cancer, by = "uniprot_symbol")

sum(str_detect(human_master$Role.in.Cancer, 'oncogene'), na.rm = TRUE)

diffs <- cancer %>%
  filter(!uniprot_symbol %in% human_master$uniprot_symbol)

stopifnot(FALSE)
write.table(x = human_master, file = "data/merged_tables/human/human_master.tsv",
            sep = "\t")