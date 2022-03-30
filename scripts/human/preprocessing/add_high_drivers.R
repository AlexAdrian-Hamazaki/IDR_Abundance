library(tidyselect)

master <- read.delim(file = 'data/merged_tables/human/actual_master.tsv', sep = "\t")

phase <- read.delim(file = 'data/phase_sep/meta/human/high_conf_drivers.txt', header = FALSE)

phase$V1 %in% master$uniprot_id

master <- master %>%
  mutate(ps_high_conf = master$uniprot_id %in% phase$V1)

typeof(master$ps_high_conf)

binarize <- function(bool) {
  if (bool == TRUE) {
    return (1)
  } else {
    return (NA)
  }
}
master$ps_high_conf <- unlist(lapply(master$ps_high_conf, binarize))

typeof(master$ps_high_conf)
master$ps_high_conf <-  as.double(as.integer(master$ps_high_conf))

write(x = master, file = "data/merged_tables/human/actual_master.tsv", sep = "\t")

saveRDS(master,file = "data/merged_tables/human/actual_master.rds" )
gg <- readRDS(file = 'data/merged_tables/human/actual_master.rds')