library(tidyverse)
big <- read.delim(file = "data/phase_sep/drlps/all_phase_sep.tsv")

yeast_scaffolds <- big %>%
  filter(Species == 'Saccharomyces cerevisiae') %>%
  filter(LLPS.Type == "Scaffold")

yeast_phase <- big %>%
  filter(Species == 'Saccharomyces cerevisiae')

yeast_clients <- big %>%
  filter(Species == 'Saccharomyces cerevisiae') %>%
  filter(LLPS.Type == "Client")

write_delim(yeast_scaffolds, file = 'data/phase_sep/drlps/yeast_drivers.tsv')

master <- readRDS(file = 'data/merged_tables/yeast/actual_master.rds')

drllps_drivers <- yeast_scaffolds %>%
  rename('uniprot_id' = "UniProt.ID") %>%
  select(uniprot_id)

drllps_drivers$uniprot_id %in% master$uniprot_id

master <- master %>%
  mutate(drlps_drivers = master$uniprot_id %in% drllps_drivers$uniprot_id)
master <- master %>%
  mutate(drpls_clients = master$uniprot_id %in% yeast_clients$UniProt.ID)
master <- master %>%
  mutate(dr_all_phase_sep = master$uniprot_id %in% yeast_phase$UniProt.ID)

binarize <- function(bool) {
  if (bool == TRUE) {
    return (1)
  } else {
    return (NA)
  }
}



master$drlps_drivers <- as.integer(unlist(lapply(master$drlps_drivers, binarize)))
master$dr_all_phase_sep <- as.integer(unlist(lapply(master$dr_all_phase_sep, binarize)))
master$drpls_clients <- as.integer(unlist(lapply(master$drpls_clients, binarize)))


write_delim( x = master, file = "data/merged_tables/yeast/actual_master.tsv", delim = "\t")

