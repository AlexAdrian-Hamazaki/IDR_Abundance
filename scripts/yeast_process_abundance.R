library(tidyr)
library(stringr)
library(dplyr)

# Purpose:
# Add the protein symbols to the abundance tables from paxDB for yeast

# Full uniprot 2 paxdb  mappings downloaded from
#https://pax-db.org/download (its the Uniprot mappings one)


# From full_uniprot_2_paxdb.tsv:
# Yeast proteins isoldated with: grep -nhr "YEASX" full* > yeast_mappings.txt
# Human proteins isolated with : grep -nhr "HUMAN" full* > human_mappings.txt


# Yeast mapping
# Reading in abundance and mapping data frames
yea_abun <- read.delim(file = "data/abundance/paxdb_abundances_integrated_whole_yeast.tsv",
                       sep = "\t",
                       header = TRUE)
yea_mappings <- read.delim(file = "data/abundance/yeast_mappings.txt",
                       header = FALSE)
colnames(yea_mappings) <- c("1", "symbol", "string_external_id", "4", '5')

yea_abun$string_external_id <- str_extract(yea_abun$string_external_id, pattern = "(?<=\\.).*")

#merge yea_mappings onto yea_abun because yea_abun is larger
yea_merged <- left_join(yea_abun, yea_mappings, by = "string_external_id", keep = FALSE, copy = FALSE)

#there's lots of duplicates, and there's lots of "proteins" in yea_abun that didn't map to anything

yea_with_symbol <- filter(yea_merged, !is.na(symbol))

print(paste(nrow(yea_with_symbol), "proteins have abundance scores"))

sum(duplicated(yea_with_symbol))

# Remove weird stuff before and after the protein symbol
yea_with_symbol <- select(yea_with_symbol, c(symbol, abundance))
symbol <- str_split(yea_with_symbol$symbol, pattern = "\\|")
symbol <- lapply(symbol, function(x) str_extract(x, pattern = ".*(?=_)"))
symbol <- lapply(symbol, function(x) return(x[2]))
yea_with_symbol$symbol <- unlist(symbol)
#AMYH
write.table(x = yea_with_symbol,
            file = "data/abundance/yea_with_symbol.tsv",
            sep = "\t",
            col.names = TRUE,
            quote = FALSE)
