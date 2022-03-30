#Purpose:

# The HI data is encoded. This script will create a HI_bin which we can merge onto human_total via uniprot_id
# Encoding is as follows
# 3,2,1 = Encoded as 'HI'
# 30, 40 = Encoded as "nHI"
# 0 or NA = Encoded as "missing"

#This script also generates summary statistics for the HI data

#Load the human_total data
human_total <- read.delim(file = 'data/merged_tables/human/human_total_data.tsv')

# split, do, merge
human_hits <- human_total %>%
  filter(HI_score == 1 | HI_score == 2 |HI_score == 3) %>%
  mutate(HI_bin = "HI")
human_fails <- human_total %>%
  filter(HI_score == 30 | HI_score == 40) %>%
  mutate(HI_bin = "nHI")
human_missing <- human_total %>%
  filter(is.na(HI_score) | HI_score == 0)%>%
  mutate(HI_bin = "missing")
stopifnot(nrow(human_hits) + nrow(human_fails) + nrow(human_missing) == nrow(human_total))

human_total_rebind <- rbind(human_hits, human_fails, human_missing)%>%
  select(uniprot_id, HI_bin)

#Write human_total_rebind. This will just be something to add onto human totaldata.tsv
write.table(x = human_total_rebind,
            file = "data/merged_tables/human/human_HI_bins.tsv",
            sep = "\t")

# See diff HI descs
HI_metadata <- human_total %>%
  select(c('HI_score', 'HI_desc')) %>%
  filter(!is.na(HI_score)) %>%
  distinct(HI_score, .keep_all = TRUE)

write.table(HI_metadata,
            file = 'data/dose/HI_metadata.tsv',
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

# Also get some summary stats for HI
HI_summary <- human_total %>%
  select(HI_desc, HI_score)%>%
  filter(!is.na(HI_desc))

no_evidence <- sum(HI_summary$HI_score == 0)
little_evidence <- sum(HI_summary$HI_score ==1  )
some_evidence <- sum(HI_summary$HI_score ==2 )
sufficient_evidence <- sum(HI_summary$HI_score ==3 )
autosomal_recessive<- sum(HI_summary$HI_score ==30 )
unlikely <- sum(HI_summary$HI_score ==40 )


HI_stats <- data.frame(names = factor(c('no_evidence',
                                     'little_evidence',
                                     'some_evidence',
                                     'sufficient_evidence',
                                     'autosomal_recessive',
                                     'unlikely',
                                     'total')),
                       levels = c('no_evidence',
                                     'little_evidence',
                                     'some_evidence',
                                     'sufficient_evidence',
                                     'autosomal_recessive',
                                     'unlikely',
                                     'total'),
                       count = c(no_evidence,
                                 little_evidence,
                                 some_evidence,
                                 sufficient_evidence,
                                 autosomal_recessive,
                                 unlikely,
                                 nrow(HI_summary)))
ggplot(HI_stats) +
  geom_col(mapping = aes(x = names, y =count))