#Purpose:
#First: Merge abundance information with protein sequence. Saves in merged_tables directory
#The reason this is a struggle is because the proteome has uniprot protein codes and uniprot protein symbols; whereas
#The abundance data has strung protein codes and string protein symbols

#Then, merge on the disorder data and happloinsufficiency data

library(tidyr)
library(dplyr)


# Load relevant tables
hum_abundance <- read.delim(file = "data/abundance/hum_with_symbol.tsv")
hum_proteome <- read.delim(file = "data/proteomes/df_human_proteome.tsv")
hum_linker <- read.delim(file = 'data/abundance/human_linker.tsv')
hum_dose <- read.delim(file= 'data/dose/ClinGen_gene_curation_list_GRCh38.tsv')
hum_iupred <- read.delim(file = "data/disorder/human_iupred.tsv",
                           header = FALSE)
hum_disopred <- read.delim(file = 'data/disorder/human_disopred3.tsv',
                             header = FALSE)

#Merge hum_abundance onto hum_linker
hum_abundance<- dplyr::rename(hum_abundance, "string_id_alias" = "symbol")
hum_merged <- left_join(df_pro_linker, y = hum_with_symbol, by = "string_id_alias")

#Merge df_human_proteome (sequence data) onto our hum_merged via uniprot_alias
df_human_proteome <- dplyr::rename(df_human_proteome, "uniprot_alias" = "symbol")
hum_merged <- left_join(hum_merged, y = df_human_proteome, by = "uniprot_alias")

#Merge disorder scores onto hum_merged via uniprot id
colnames(human_disopred) <- c("uniprot_id",'disopred')
colnames(human_iupred) <- c("uniprot_id",'iupred')
hum_merged <- left_join(hum_merged, y = human_disopred, by = 'uniprot_id')
hum_merged <- left_join(hum_merged, y = human_iupred, by = 'uniprot_id')


hum_no_abundance <- hum_merged %>%
  filter(is.na(hum_merged$abundance))
hum_no_seq<- hum_merged %>%
  filter(is.na(hum_merged$seq))

hum_abundance <- hum_merged %>%
  filter(!is.na(hum_merged$abundance) &!is.na(hum_merged$seq))

print(paste(nrow(hum_no_abundance), 'proteins in the human proteome did not have abundance values'))
print(paste(nrow(hum_no_seq), 'proteins in the human proteome did not have sequences values'))
print(paste(nrow(hum_abundance), 'proteins in the human proteome have abundance values'))

write.table(x = hum_abundance,
            file = 'data/merged_tables/human_abundanceXsequencesXdisorder.tsv',
            col.names = TRUE,
            quote = FALSE,
            sep = "\t"
)
