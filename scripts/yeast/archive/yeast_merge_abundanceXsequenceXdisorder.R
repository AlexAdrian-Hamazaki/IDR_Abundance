#Purpose:
#Merge yeast abundance information with protein sequence. Saves in merged_tables directory

library(tidyr)
library(dplyr)



yea_with_symbol <- read.delim(file = "data/abundance/yea_with_symbol.tsv")
df_yeast_proteome <- read.delim(file = "data/proteomes/df_yeast_proteome.tsv")
yeast_iupred <- read.delim(file = "data/disorder/yeast_iupred.tsv",
                           header = FALSE)
yeast_disopred <- read.delim(file = 'data/disorder/yeast_disopred3.tsv',
                             header = FALSE)

#Merge yea_with_symbols (our abundance data) onto df_pro_linker via string_id_alias
yea_with_symbol <- dplyr::rename(yea_with_symbol, "string_id_alias" = "symbol")
yea_merged <- left_join(df_pro_linker, y = yea_with_symbol, by = "string_id_alias")

#Merge df_yeast_proteome (sequence data) onto our yea_merged via uniprot_alias
df_yeast_proteome <- dplyr::rename(df_yeast_proteome, "uniprot_alias" = "symbol")
yea_merged <- left_join(yea_merged, y = df_yeast_proteome, by = "uniprot_alias")

#Merge disorder scores onto yea_merged via uniprot id
colnames(yeast_disopred) <- c("uniprot_id",'disopred')
colnames(yeast_iupred) <- c("uniprot_id",'iupred')
yea_merged <- left_join(yea_merged, y = yeast_disopred, by = 'uniprot_id')
yea_merged <- left_join(yea_merged, y = yeast_iupred, by = 'uniprot_id')


yea_no_abundance <- yea_merged %>%
  filter(is.na(yea_merged$abundance))
yea_no_seq<- yea_merged %>%
  filter(is.na(yea_merged$seq))

yea_abundance <- yea_merged %>%
  filter(!is.na(yea_merged$abundance) &!is.na(yea_merged$seq))

print(paste(nrow(yea_no_abundance), 'proteins in the yeast proteome did not have abundance values'))
print(paste(nrow(yea_no_seq), 'proteins in the yeast proteome did not have sequences values'))
print(paste(nrow(yea_abundance), 'proteins in the yeast proteome have abundance values'))

write.table(x = yea_abundance,
            file = 'data/merged_tables/yeast_abundanceXsequencesXdisorder.tsv',
            col.names = TRUE,
            quote = FALSE,
            sep = "\t"
)
