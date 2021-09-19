library(dplyr)
library(stringr)


dose <- read.delim(file = "data/dose/processed_doses.tsv")
yea_proteins <-  read.delim(file = 'data/merged_tables/yeast_abundanceXsequencesXdisorder.tsv')

yea_dose_proteins_uni <- left_join(x = yea_proteins, y = dose, by = c('uniprot_alias' = "symbol"))
yea_dose_proteins_uni <- yea_dose_proteins_uni%>%
  filter(!is.na(HI_score))
yea_dose_proteins_string <- left_join(x = yea_proteins, y = dose, by = c('string_id_alias' = 'symbol'))
yea_dose_proteins_string <- yea_dose_proteins_string%>%
  filter(!is.na(HI_score))

dose_sensitive <- rbind(yea_dose_proteins_uni, yea_dose_proteins_string)
dose_sensitive <- dose_sensitive %>%
  unique(by = "uniprot_id")