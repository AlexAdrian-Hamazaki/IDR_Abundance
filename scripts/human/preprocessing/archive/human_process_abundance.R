library(tidyr)
library(stringr)
library(dplyr)

# Purpose:
# Add the protein symbols to the abundance tables from paxDB for human

#See data/abundance/download_links.txt for the URLs of where the abundance and mapping data was downloaded from

# human mapping
# Reading in abundance and mapping data frames
hum_abun <- read.delim(file = "../../../../data/abundance/human/paxdb_abundances_integrated_whole_human.tsv",
                       sep = "\t",
                       header = TRUE)
hum_mappings <- read.delim(file = "../../../data/abundance/9606.protein.info.v11.5.txt",
                           header = TRUE,
                           sep = "\t",
                           quote = "",
                           fill = FALSE)

hum_mappings <- rename(hum_mappings, "string_external_id" = "X.string_protein_id")
hum_mappings <- rename(hum_mappings, "symbol" = "preferred_name")


expected_complete <- sum(hum_abun$string_external_id %in% hum_mappings$string_external_id)
print(paste(expected_complete, 'proteins with abundance scores are expected to have protein symbols'))


#merge hum_mappings onto hum_abun because hum_abun is larger
hum_merged <- left_join(hum_abun, hum_mappings, by = "string_external_id", keep = FALSE)

#there's lots of duplicates, and there's lots of "proteins" in hum_abun that didn't map to anything

hum_with_symbol <- filter(hum_merged, !is.na(symbol))
hum_without_symbol <- filter(hum_merged, is.na(symbol))
print(paste(nrow(hum_without_symbol), "proteins with abundance scores do not have symbols"))
print(paste(nrow(hum_with_symbol), "proteins with abundance scores have symbols"))
print(paste(expected_complete/nrow(hum_abun), '% of human genes have abundance scores'))

# Hard to tell why there arn't symbols for a lot of the proteins with abundance - But I guess we just don't have abundance data for
# all human genes

########## Looking at duplicates
# sum(duplicated(hum_with_symbol$symbol))
# unique(hum_with_symbol$symbol)
# hum_with_symbol[(hum_with_symbol$symbol),]
#There are no duplicate symbols,but there are "NA's" that are duplicated


# Remove weird stuff before and after the protein symbol
hum_with_symbol <- select(hum_with_symbol, c(string_external_id, symbol, abundance, protein_size))

#Save
write.table(x = hum_with_symbol,
            file = "../../../../data/abundance/hum_with_symbol.tsv",
            sep = "\t",
            col.names = TRUE,
            quote = FALSE)
