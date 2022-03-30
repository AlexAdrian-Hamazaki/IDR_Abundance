library(dplyr)
library(tidyr)

dose <- read.delim(file = "data/dose/ClinGen_gene_curation_list_GRCh38.tsv")

#Remove genes that are autosomal recessive(30) or do not have evidence(0) and have little evidence(1) of dose sensitive pathogenicity

pro_dose <- dose %>%
  filter(Haploinsufficiency.Score != 0)%>%
  dplyr::select(c("Gene.Symbol","Haploinsufficiency.Score"))
colnames(pro_dose) <- c('symbol', 'HI_score')

write.table(x = pro_dose,
            file = "data/dose/processed_doses.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)