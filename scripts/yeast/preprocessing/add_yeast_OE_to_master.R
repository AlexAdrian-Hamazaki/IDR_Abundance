library(tidyselect)

master <- read.delim(file = 'data/merged_tables/yeast/actual_master.tsv',
                     sep = '\t')
ff <- data.frame(ff)
typeof(ff)
ff$new_oe <- as.integer(ff$new_oe)
typeof(ff$new_oe)
typeof(ff$HI)

write(ff, file = "data/merged_tables/yeast/actual_master.tsv", sep = "\t")


new_OE <- read.delim(file = 'data/dose/yeast/new_OE.csv', sep = "\t")
new_OE$uniprot_id
new_OE$uniprot_id %in% master$uniprot_id

master2 <- master %>%
  mutate(new_oe = master$uniprot_id %in% new_OE$uniprot_id)

acc <- list()

binarize <- function(bool) {
  if (bool == TRUE) {
    return (1)
  } else {
    return (NA)
  }
}
ff <- lapply(master2$new_oe, binarize)
unlist(ff)
master2$new_oe <- as.integer(unlist(ff))

write(master2, file = "data/merged_tables/yeast/actual_master.tsv",
      sep = "\t")

saveRDS(master2, file = "data/merged_tables/yeast/actual_master.rds")
gg <- readRDS('data/merged_tables/yeast/actual_master.rds')