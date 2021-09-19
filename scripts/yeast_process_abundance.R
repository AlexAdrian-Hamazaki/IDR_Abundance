library(tidyr)
library(stringr)
library(dplyr)

# Purpose:
# Add the protein symbols to the abundance tables from paxDB for yeast

#See data/abundance/download_links.txt for the URLs of where the abundance and mapping data was downloaded from

# Yeast mapping
# Reading in abundance and mapping data frames
yea_abun <- read.delim(file = "data/abundance/paxdb_abundances_integrated_whole_yeast.tsv",
                       sep = "\t",
                       header = TRUE)
yea_mappings <- read.delim(file = "data/abundance/4932.protein.info.v11.5.txt",
                       header = TRUE)
yea_mappings <- rename(yea_mappings, "string_external_id" = "X.string_protein_id")
yea_mappings <- rename(yea_mappings, "symbol" = "preferred_name")


expected_complete <- sum(yea_abun$string_external_id %in% yea_mappings$string_external_id)
print(paste(expected_complete, 'proteins with abundance scores are expected to have protein symbols"'))


#merge yea_mappings onto yea_abun because yea_abun is larger
yea_merged <- left_join(yea_abun, yea_mappings, by = "string_external_id", keep = FALSE)

#there's lots of duplicates, and there's lots of "proteins" in yea_abun that didn't map to anything

yea_with_symbol <- filter(yea_merged, !is.na(symbol))
yea_without_symbol <- filter(yea_merged, is.na(symbol))
print(paste(nrow(yea_without_symbol), "proteins with abundance scores do not have symbols"))
print(paste(nrow(yea_with_symbol), "proteins with abundance scores have symbols"))

########## Looking at duplicates
# sum(duplicated(yea_with_symbol$symbol))
# unique(yea_with_symbol$symbol)
# yea_with_symbol[(yea_with_symbol$symbol),]
#There are no duplicate symbols,but there are "NA's" that are duplicated


# Remove weird stuff before and after the protein symbol
yea_with_symbol <- select(yea_with_symbol, c(symbol, abundance, protein_size))

#Save
write.table(x = yea_with_symbol,
            file = "data/abundance/yea_with_symbol.tsv",
            sep = "\t",
            col.names = TRUE,
            quote = FALSE)
