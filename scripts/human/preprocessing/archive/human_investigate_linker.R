library(dplyr)
library(stringr)

linker <- read.delim("data/abundance/human_raw_linker.tsv")

linker_9606 <- linker %>%
  filter(str_detect(uniprot_alias, pattern = "9606"))