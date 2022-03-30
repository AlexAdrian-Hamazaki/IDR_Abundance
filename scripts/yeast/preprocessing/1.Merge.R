#Purpose:

# Get a merged table that merges
# 1. Abundance info. 2. HI info. 3. OE info and 4. Disorder info
library(tidyr)
library(dplyr)

#~~~~~~~~~ Processing Disorder Info ~~~~~~~~~~~~
disopred <- read.delim(file = "data/disorder/yeast/yeast_disopred3.tsv", header = FALSE)
iupred <- read.delim(file = "data/disorder/yeast/yeast_iupred.tsv", header = FALSE)

colnames(disopred) <- c("uniprot_id", "disopred3perc")
colnames(iupred) <- c("uniprot_id", "iupredperc")


#~~~~~~~~~ Processing Abundance Info ~~~~~~~~~~~~
abundance <- read.delim(file = "data/abundance/yeast/paxdb_abundances_integrated_whole_yeast.tsv",
                        sep = "\t",
                        header = TRUE)
strdb_to_uniprot <- read.delim(file = "data/abundance/yeast/string_to_uniprot.txt")

abundance_pro <- left_join(abundance, strdb_to_uniprot, by = c("string_external_id" = "yourlist.M2021110292C7BAECDB1C5C413EE0E0348724B6822776773"))

message(paste("Removed", sum(is.na(abundance_pro$Entry)), " abundance genes because we could not translate string ID to uniprot"))

abundance_pro <- abundance_pro %>%
  filter(!is.na(Entry)) %>%
  select(Entry, abundance, "Gene.names") %>%
  rename("uniprot_symbols" = "Gene.names") %>%
  rename("uniprot_id" = "Entry")



#~~~~~~~~~ Processing HI Info ~~~~~~~~~~~~
HI <- read.delim(file = "data/dose/yeast/yeast_HI_genes.tsv",
                 sep = "\t")
HI <- HI %>%
  select("ORF") %>%
  rename("entrez_id" = "ORF") %>%
  mutate(HI = 1)

#~~~~~~~~~ Processing OE ~~~~~~~~~~~~
OE <- read.delim(file = "data/dose/yeast/yeast_OE_Genes.txt",
                   sep = "\t")
OE <- OE %>%
  select("Systematic.Name") %>%
  rename("entrez_id" = "Systematic.Name") %>%
  mutate(OE = 1)

#~~~~~~~~~ Process entrez_to_uniprot.tsv for merging~~~~~~~~~~~~
entrez_to_uniprot <- read.delim(file = "data/dose/yeast/entrez_to_uniprot.tsv",
                                sep = "\t")
colnames(entrez_to_uniprot) <- c("entrez_id", "uniprot_id")
entrez_to_uniprot <- entrez_to_uniprot[,1:2]
#~~~~~~~~~ Building Merged Data ~~~~~~~~~~~~


merge <- disopred
merge <- left_join(merge, iupred, by = "uniprot_id")
merge <- left_join(merge, abundance_pro, by = 'uniprot_id')
merge <- left_join(merge, entrez_to_uniprot, by = "uniprot_id")
merge <- left_join(merge, HI, by = "entrez_id")
merge <- left_join(merge, OE, by = "entrez_id")

dir.create("data/merged_tables/yeast", showWarnings = FALSE)
write.table(merge, file = "data/merged_tables/yeast/actual_master.tsv",
            sep = "\t")
