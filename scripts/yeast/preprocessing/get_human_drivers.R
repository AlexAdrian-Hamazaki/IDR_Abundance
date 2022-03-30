library(tidyverse)
big <- read.delim(file = "data/phase_sep/drlps/all_phase_sep.tsv")

human_scaffolds <- big %>%
  filter(Species == 'Homo sapiens') %>%
  filter(LLPS.Type == "Scaffold")

human_phase <- big %>%
  filter(Species == 'Homo sapiens')

write_delim(human_scaffolds, file = 'data/phase_sep/drlps/human_drivers.tsv')

master <- readRDS(file = 'data/merged_tables/human/actual_master.rds')

drllps_drivers <- human_scaffolds %>%
  rename('uniprot_id' = "UniProt.ID") %>%
  select(uniprot_id)

drllps_drivers$uniprot_id %in% master$uniprot_id

master <- master %>%
  mutate(drlps_drivers = master$uniprot_id %in% drllps_drivers$uniprot_id)

binarize <- function(bool) {
  if (bool == TRUE) {
    return (1)
  } else {
    return (NA)
  }
}
master <- master %>%
  mutate(dr_all_phase_sep = master$uniprot_id %in% human_phase$UniProt.ID)

clients <- human_phase %>%
  filter(LLPS.Type == "Client")

master$drpls_clients <- master$uniprot_id %in% clients$UniProt.ID


master$drlps_drivers <- as.integer(unlist(lapply(master$drlps_drivers, binarize)))
master$dr_all_phase_sep <- as.integer(unlist(lapply(master$dr_all_phase_sep, binarize)))
master$drpls_clients <- as.integer(unlist(lapply(master$drpls_clients, binarize)))



write_delim( x = master, file = "data/merged_tables/human/actual_master.tsv", delim = "\t")

