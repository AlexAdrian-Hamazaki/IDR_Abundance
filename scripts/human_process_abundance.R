library(tidyr)
library(stringr)
library(dplyr)

# Purpose:
# Add the protein symbols to the abundance tables from paxDB for human

# Full uniprot 2 paxdb  mappings downloaded from
#https://pax-db.org/download (its the Uniprot mappings one)


# From full_uniprot_2_paxdb.tsv:
# Human proteins isolated with : grep -nhr "HUMAN" full* > human_mappings.txt


# human mapping
# Reading in abundance and mapping data frames
hum_abun <- read.delim(file = "data/abundance/paxdb_abundances_integrated_whole_human.tsv",
                       sep = "\t",
                       header = TRUE)
hum_mappings <- read.delim(file = "data/abundance/human_mappings.txt",
                       header = FALSE)
colnames(hum_mappings) <- c("1", "symbol", "string_external_id", "4", '5')

hum_abun$string_external_id <- str_extract(hum_abun$string_external_id, pattern = "(?<=\\.).*")

#merge hum_mappings onto hum_abun because hum_abun is larger
hum_merged <- left_join(hum_abun, hum_mappings, by = "string_external_id", keep = FALSE, copy = FALSE)

#there's lots of duplicates, and there's lots of "proteins" in hum_abun that didn't map to anything

hum_with_symbol <- filter(hum_merged, !is.na(symbol))

print(paste(nrow(hum_with_symbol), "proteins have abundance scores"))

sum(duplicated(hum_with_symbol))

# Remove weird stuff before and after the protein symbol
hum_with_symbol <- select(hum_with_symbol, c(symbol, abundance))
symbol <- str_split(hum_with_symbol$symbol, pattern = "\\|")
symbol <- lapply(symbol, function(x) str_extract(x, pattern = ".*(?=_)"))
symbol <- lapply(symbol, function(x) return(x[2]))
hum_with_symbol$symbol <- unlist(symbol)


write.table(x = hum_with_symbol,
            file = "data/abundance/hum_with_symbol.tsv",
            sep = "\t",
            col.names = TRUE,
            quote = FALSE)
