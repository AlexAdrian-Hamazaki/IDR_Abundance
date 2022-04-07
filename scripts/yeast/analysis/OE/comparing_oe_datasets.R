library(tidyverse)

master <- read.delim(file = 'data/merged_tables/yeast/actual_master.tsv', sep = "\t")

oe <- master %>%
  select(uniprot_id, OE, new_oe)

both <- oe%>%
  filter(!is.na(OE)) %>%
  filter(!is.na(new_oe))

sum(oe$OE, na.rm =  TRUE)