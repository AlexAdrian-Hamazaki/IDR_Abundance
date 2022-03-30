#Purpose:

# Get a merged table that merges 1. Abundance info. 2. HI info. 3. OE info and 4. Disorder info
library(tidyr)
library(dplyr)
library(stringr)

#~~~~~~~~~ Processing Disorder Info ~~~~~~~~~~~~
human_disopred <- read.delim(file = "data/disorder/human/human_disopred3.tsv", header = FALSE)
human_iupred <- read.delim(file = "data/disorder/human/human_iupred.tsv", header = FALSE)

colnames(human_disopred) <- c("uniprot_id", "disopred3perc")
colnames(human_iupred) <- c("uniprot_id", "iupredperc")


#~~~~~~~~~ Processing Abundance Info ~~~~~~~~~~~~
abundance <- read.delim(file = "data/abundance/human/paxdb_abundances_integrated_whole_human.tsv",
                        sep = "\t",
                        header = TRUE)
strdb_to_uniprot <- read.delim(file = "data/abundance/human/string_to_uniprot.txt")

abundance_pro <- left_join(abundance, strdb_to_uniprot, by = c("string_external_id" = "yourlist.M20211025A084FC58F6BBA219896F365D15F2EB4424DADCH"))

message(paste("Removed", sum(is.na(abundance_pro$Entry)), " abundance genes because we could not translate string ID to uniprot"))
abundance_pro <- abundance_pro %>%
  filter(!is.na(Entry)) %>%
  select(Entry, abundance, "Gene.names") %>%
  rename("uniprot_symbols" = "Gene.names") %>%
  rename("uniprot_id" = "Entry")



#~~~~~~~~~ Processing HI Info ~~~~~~~~~~~~
HI <- read.delim(file = "data/dose/human/ClinGen_gene_curation_list_GRCh38.tsv",
                 sep = "\t")
HI <- HI %>%
  select("Gene.ID", "Haploinsufficiency.Score") %>%
  rename("entrez_id" = "Gene.ID") %>%
  rename("HI" = "Haploinsufficiency.Score")
HI_true <- HI %>%
  filter(HI %in% c(1,2,3))%>%
  mutate(is_HI = 1)
HI_false <- HI %>%
  filter(!HI %in% c(1,2,3))%>%
  mutate(is_HI = 0)
HI <- rbind(HI_true, HI_false) %>%
  select(entrez_id, is_HI)

HI$entrez_id <- as.character(HI$entrez_id)

#~~~~~~~~~ Processing oncogene and HI Info ~~~~~~~~~~~~
Onco <- read.delim(file = "data/dose/human/Census_allMon Oct 25 22 33 06 2021.tsv",
                   sep = "\t")
Onco <- Onco %>%
  select("Entrez.GeneId", "Role.in.Cancer") %>%
  rename("entrez_id" = "Entrez.GeneId") %>%
  rename("cancer_role" = "Role.in.Cancer")
Onco$entrez_id <- as.character(Onco$entrez_id)

Onco_true <- Onco %>%
  filter(str_detect(string = cancer_role, pattern = "oncogene")) %>%
  mutate(oncogene = 1)
Onco_false <- Onco %>%
  filter(!str_detect(string = cancer_role, pattern = "oncogene")) %>%
  mutate(oncogene = 0)
Onco <- rbind(Onco_true, Onco_false) %>%
  select('entrez_id', 'oncogene')


#~~~~~~~~~ Process entrez_to_uniprot.tsv for merging~~~~~~~~~~~~
entrez_to_uniprot <- read.delim(file = "data/dose/human/entrez_to_uniprot.tsv",
                                sep = "\t")
colnames(entrez_to_uniprot) <- c("entrez_id", "uniprot_id")
#~~~~~~~~~ Building Merged Data ~~~~~~~~~~~~


merge <- human_disopred
merge <- left_join(merge, human_iupred, by = "uniprot_id")
merge <- left_join(merge, abundance_pro, by = 'uniprot_id')
merge <- left_join(merge, entrez_to_uniprot, by = "uniprot_id")
merge <- left_join(merge, HI, by = "entrez_id")
merge <- left_join(merge, Onco, by = "entrez_id")

write.table(merge, file = "data/merged_tables/human/actual_master.tsv",
            sep = "\t")

stopifnot(FALSE)

#~~~~~~~~~ Investigating Oncogene Merging ~~~~~~~~~~~~
#Gotta see how many oncogenes actually mapped.
sum(!is.na(merge$cancer_role))
# Duplicates
cancer_test <- filter(merge, !is.na(cancer_role))
cancer_test[duplicated(cancer_test$entrez_id),]

#Note: entrez_ids 2778, 1523, and 7307, 1029 all map to more than 1 uniprot_ids. So when merging the oncogenes onto the merge
# table, you will actually have 'more' oncogenes than you can possibly merge (due to duplications). However,
# for each of the entrez ids, only 1 of each actually have abundance information; thus, the duplicates will be removed
# and we won't have false positive oncogene data.

# There is also a problem merging 3020: #TODO
