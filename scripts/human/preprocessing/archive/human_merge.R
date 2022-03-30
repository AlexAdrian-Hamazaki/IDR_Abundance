#Purpose:
#First: Merge abundance information with protein sequence. Saves in merged_tables directory
#The reason this is a struggle is because the proteome has uniprot protein codes and uniprot protein symbols; whereas
#The abundance data has strung protein codes and string protein symbols

#Then, merge on the disorder data and happloinsufficiency data

library(tidyr)
library(dplyr)


# Load relevant tables
hum_abundance <- read.delim(file = "../../../../data/abundance/hum_with_symbol.tsv")
hum_dose <- read.delim(file= '../../../data/dose/ClinGen_gene_curation_list_GRCh38.tsv')
hum_iupred <- read.delim(file = "../../../data/disorder/human_iupred.tsv",
                           header = FALSE)
hum_disopred <- read.delim(file = '../../../data/disorder/human_disopred3.tsv',
                             header = FALSE)
hum_linker <- read.delim(file = '../../../../data/linkers/human_pro_linker.tsv')

# We will be merging ALL of our different information onto hum_linker because it has all of the possible codes and ids we need
# We will be merging onto final_merge

# Merge hum_abundance onto hum_linker via string_id
hum_abundance <- hum_abundance %>%
  rename("string_id" = "string_external_id") %>%
  select(c('string_id', 'abundance', 'symbol'))

final_merge <- left_join(hum_linker, y = hum_abundance, by = "string_id")
stopifnot(all(final_merge$stringdb_symbol == final_merge$symbol, na.rm= TRUE))
# If programms fails this checkpoint then there is a huge problem with your strindb symbols that needs investigating

final_merge <- select(final_merge, c(-symbol))

# Merge hum_dose onto linker's string db and uniprot db symbols. then add that dose_merge data frame onto final_merge)'
hum_dose <- hum_dose %>%
  rename('dose_symbol' = 'Gene.Symbol') %>%
  rename('HI_score' = 'Haploinsufficiency.Score') %>%
  rename("HI_desc" = 'Haploinsufficiency.Description') %>%
  rename('TRIPLO_score' = "Triplosensitivity.Score") %>%
  rename('TRIPLO_desc' = 'Triplosensitivity.Description') %>%
  rename('OMIM_ID' = "Loss.phenotype.OMIM.ID") %>%
  select('dose_symbol',
         'HI_score',
         'HI_desc',
          'TRIPLO_score',
         'TRIPLO_desc',
         'OMIM_ID')
dose_merge_uni <- left_join(hum_linker, hum_dose, by = c('uniprot_symbol' = 'dose_symbol'), keep = TRUE)
dose_merge_uni <- dose_merge_uni %>%
  select('uniprot_id',
         'HI_score',
         'HI_desc',
          'TRIPLO_score',
         'TRIPLO_desc',
         'OMIM_ID') %>%
  filter(!is.na(HI_score))

dose_merge_string <- left_join(hum_linker, hum_dose, by = c('stringdb_symbol' = 'dose_symbol'))
dose_merge_string <- dose_merge_string %>%
  select('uniprot_id',
         'HI_score',
         'HI_desc',
          'TRIPLO_score',
         'TRIPLO_desc',
         'OMIM_ID') %>%
  filter(!is.na(HI_score))

dose_merge_combined <- rbind(dose_merge_uni, dose_merge_string) %>%
  distinct(uniprot_id, .keep_all = TRUE)

final_merge <- left_join(final_merge, dose_merge_combined, by = "uniprot_id")

message(paste(sum(!is.na(final_merge$HI_score)), ": Dose sensiive scores were able to map"))

#Merge disorder scores onto hum_merged via uniprot id
colnames(hum_disopred) <- c("uniprot_id",'disopred')
colnames(hum_iupred) <- c("uniprot_id",'iupred')
final_merge <- left_join(final_merge, y = hum_disopred, by = 'uniprot_id')
final_merge <- left_join(final_merge, y = hum_iupred, by = 'uniprot_id')

#Remove some duplicate Uniprot IDs
final_merge %>% distinct(uniprot_symbol)

final_merge <- final_merge[,c(1,2,3,4,7,13,14,8,9,10,11,12,5,6)]

prep_matrix <- list(as.numeric(final_merge[,'abundance']) ,
                    as.numeric(final_merge[, 'disopred']),
                    as.numeric(final_merge[, 'iupred']),
                    as.numeric(final_merge[, 'HI_score']),
                    as.numeric(final_merge[, 'TRIPLO_score'])
)


final_matrix <- matrix(do.call(cbind, prep_matrix),
                       ncol = 5,
                       nrow = nrow(final_merge),
                       dimnames = list(c(final_merge$uniprot_symbol),c('abundance',
                                             'dispored',
                                             'iupred',
                                             'HI_score',
                                              'TRIPLO_score')))
final_metadata <- final_merge %>%
  select('uniprot_id',
         'uniprot_symbol',
         'string_id',
         'stringdb_symbol',
         'HI_desc',
         'TRIPLO_desc',
         'OMIM_ID',
         'protein_size',
         'seq')
rownames(final_metadata) <- final_metadata$uniprot_symbol


write.table(x = final_merge,
            file = '../../../../data/merged_tables/human/human_total_data.tsv',
            col.names = TRUE,
            quote = FALSE,
            sep = "\t"
)

#Extract only matrix data

# Extract only metadata
stopifnot(FALSE)
hum_no_abundance <- hum_merged %>%
  filter(is.na(hum_merged$abundance))
hum_no_seq<- hum_merged %>%
  filter(is.na(hum_merged$seq))

hum_abundance <- hum_merged %>%
  filter(!is.na(hum_merged$abundance) &!is.na(hum_merged$seq))

print(paste(nrow(hum_no_abundance), 'proteins in the human proteome did not have abundance values'))
print(paste(nrow(hum_no_seq), 'proteins in the human proteome did not have sequences values'))
print(paste(nrow(hum_abundance), 'proteins in the human proteome have abundance values'))
