# Purpose:

# See if Highly Disordered, Low Abundant proteins are enriched for delatorious overexpression phenotypes (via oncogene phenotype)
# compared to the proteome

# looking at Disopred3% score

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

source(file = "scripts/functions.R")


#~~~~~~~~~ Open Master ~~~~~~~~~~~~

master <- read.delim('data/merged_tables/human/actual_master.tsv',
                     sep = "\t")
# log2+1 abundance
master$abundance <- log2(master$abundance+1)

#~~~~~~~~~ Global Parameters ~~~~~~~~~~~~
#Abundnce percent is whatever percent you want to select as "low" abundant. Must be in 0.4 format not 40%.
abundance_percent <- 0.4

#disorder percent is whatever percent you want to select as  "high" disorder Must be in 0.4 format not 40%.
disorder_percent <- 0.35

print(paste("Abundance Percent:", abundance_percent))
print(paste("Disorder Percent:", disorder_percent))


#~~~~~~~~~ Get all low abundant proteins in proteome ~~~~~~~~~~~~

LA <- master %>%
  filter(!is.na(abundance) & !is.na(disopred3perc)) %>%
  get_highest_percent(abundance, 0.33, from = 'low')


LA_counts <- LA %>%
    get_int_counts('oncogene')


#~~~~~~~~~ Of get High Disorder low/mid/high abundant proteins ~~~~~~~~~~~~

HD_LA <- master %>%
  filter(!is.na(disopred3perc) & !is.na(abundance)) %>%
  get_highest_percent(disopred3perc, disorder_percent) %>%
  get_highest_percent(abundance, abundance_percent, from = 'low')

HD_MA <- master %>%
  filter(!is.na(disopred3perc) & !is.na(abundance)) %>%
  get_highest_percent(disopred3perc, disorder_percent) %>%
  get_middle_percent(abundance, abundance_percent, from = 'low')

HD_HA <- master %>%
  filter(!is.na(disopred3perc) & !is.na(abundance)) %>%
  get_highest_percent(disopred3perc, disorder_percent) %>%
  get_lowest_percent(abundance, abundance_percent, from = 'low')

#~~~~~~~~~ Count number of oncogenes in each  ~~~~~~~~~~~~

HD_LA_counts <- HD_LA %>%
    get_int_counts('oncogene')

HD_MA_counts <- get_int_counts(HD_MA,'oncogene')

HD_HA_counts <- get_int_counts(HD_HA,'oncogene')


#~~~~~~~~~ Bring contingency table together and do fisher exact test ~~~~~~~~~~~~

fisher_table <- rbind(LA_counts, HD_LA_counts, HD_MA_counts, HD_HA_counts)
rownames(fisher_table) <- c("LA", "HD_LA", "HD_MA", "HD_HA")

fisher_HD_LA <- fisher_table %>%
  filter(rownames(fisher_table) %in% c("LA", "HD_LA")) %>%
  select(-c(total))%>%
  fisher.test()
print(paste("P.val: HD_LA vs LA:", fisher_HD_LA$p.value))


fisher_HD_MA <- fisher_table %>%
  filter(rownames(fisher_table) %in% c("LA", "HD_MA")) %>%
  select(-c(total))%>%
  fisher.test()
print(paste("P.val: HA_MA vs LA:", fisher_HD_MA$p.value))


fisher_HD_HA <- fisher_table %>%
  filter(rownames(fisher_table) %in% c("LA", "HD_HA")) %>%
  select(-c(total))%>%
  fisher.test()
print(paste("P.val: HD_HA vs LA:", fisher_HD_HA$p.value))

#Mereley having HA makes it significant

fisher_elongate <- fisher_table%>%
  select(-c(total)) %>%
  mutate(names  = rownames(fisher_table)) %>%
    melt()
name_levels <- c('HD_HA', 'HD_MA', 'HD_LA', "LA")
ggplot(data = fisher_elongate) +
  geom_bar(mapping =  aes(x = factor(names, levels = name_levels), y = value, fill = variable),
           position='stack', stat = 'identity')+
  labs(title = "HI counts for HD proteins across different abundance levels",
       x = "Different Abundance Levels",
       y = "counts",
      fill = "HI Protein Presence")+
  scale_fill_discrete(labels = c('HI Proteins', 'Not HI Proteins'))+
  scale_x_discrete(labels = name_levels)

#~~~~~~~~~ Comparing Entire proteoms vs High Disorder for OE enrichment ~~~~~~~~~~~~
proteome <- master %>%
  filter(!is.na(disopred3perc) & !is.na(abundance))
proteome_counts <- get_int_counts(proteome, 'oncogene')

HD <- master  %>%
  filter(!is.na(disopred3perc) & !is.na(abundance)) %>%
  get_highest_percent(disopred3perc, disorder_percent)
HD_counts <- get_int_counts(HD, 'oncogene')

vaboori <- rbind(proteome_counts, HD_counts)
rownames(vaboori) <- c('proteome', 'HD')
vaboori_fisher <- vaboori %>%
  select(-c(total)) %>%
  fisher.test()

print(paste("P.val: HD vs proteome:", vaboori_fisher$p.value))

vaboori_perc <- vaboori%>%
    summarize(percent = counts/not_counts*100,
              name = rownames(vaboori))

ggplot(vaboori_perc) +
  geom_col(mapping = aes(x = name, y = percent), fill = 'steelblue')+
  theme_classic()+
  labs(title = "Vaboori Oncogene Replication",
       y = "Percent OE",
       x = 'Proteins')

stopifnot(FALSE)

#~~~~~~~~~ Now as a control see if the high abundance high disopred is enriched with OE, I hope it isn't ~~~~~~~~~~

HD_HA <- master %>%
  filter(!is.na(disopred3perc) & !is.na(abundance)) %>%
  get_highest_percent(disopred3perc, disorder_percent) %>%
  get_highest_percent(abundance, abundance_percent)



#~~~~~~~~~ Plottong ~~~~~~~~~~~~
master_plot <- master %>%
  filter(!is.na(abundance) & !is.na(disopred3perc))
ggplot(master_plot) +
  geom_histogram(mapping = aes(abundance))

ggplot(HD_LA) +
  geom_histogram(mapping = aes(abundance))
ggplot(HD_HA) +
  geom_histogram(mapping = aes(abundance))



#~~~~~~~~~ Do fisher for lowest, mid and highest abundance ~~~~~~~~~~~~

#HD_LA is already done, now we use that percentage to bin the mid and high abundance

HD_MA <- master %>%
  filter(!is.na(disopred3perc) & !is.na(abundance)) %>%
  get_highest_percent(disopred3perc, disorder_percent) %>%
  get_middle_percent(abundance, abundance_percent, from = 'low')

HD_HA <- master %>%
  filter(!is.na(disopred3perc) & !is.na(abundance)) %>%
  get_highest_percent(disopred3perc, disorder_percent) %>%
  get_lowest_percent(abundance, abundance_percent, from = 'low')

# Get count Tables

HD_MA_counts <- HD_MA %>%
    summarize(df_rows = nrow(HD_MA),
              oncogene_counts = count_string(cancer_role, 'oncogene'),
              not_oncogene_counts = count_string(cancer_role, 'oncogene', type = 'not')
    )

HD_HA_counts <- HD_HA %>%
    summarize(df_rows = nrow(HD_HA),
              oncogene_counts = count_string(cancer_role, 'oncogene'),
              not_oncogene_counts = count_string(cancer_role, 'oncogene', type = 'not')
    )

# Bind contingency table

total_fisher_table <- rbind(LA_counts, HD_LA_counts, HD_MA_counts, HD_HA_counts)
rownames(total_fisher_table) <- c('LA', 'HD_LA', 'HD_MA', 'HD_HA')

#  Fishers for High abundance, mid abundance and low abundance

HD_HA_fisher <- total_fisher_table %>%
  select(not_oncogene_counts, oncogene_counts) %>%
  filter(rownames(total_fisher_table)%in% c('LA', 'HD_HA'))%>%
  fisher.test()
print(paste("HD_HA:", HD_HA_fisher$p.value))

HD_MA_fisher <- total_fisher_table %>%
  select(not_oncogene_counts, oncogene_counts) %>%
  filter(rownames(total_fisher_table)%in% c('LA', 'HD_MA'))%>%
  fisher.test()
print(paste("HD_MA:",HD_MA_fisher$p.value))


HD_LA_fisher <- total_fisher_table %>%
  select(not_oncogene_counts, oncogene_counts) %>%
  filter(rownames(total_fisher_table)%in% c('LA', 'HD_LA'))%>%
  fisher.test()
print(paste("HD_LA:",HD_LA_fisher$p.value))


# Do a fisher for the entirey of High Disorder

HD <- master %>%
  filter(!is.na(disopred3perc) & !is.na(abundance)) %>%
  get_highest_percent(disopred3perc, disorder_percent)

HD_counts <- HD%>%
    summarize(df_rows = nrow(HD),
              oncogene_counts = count_string(cancer_role, 'oncogene'),
              not_oncogene_counts = count_string(cancer_role, 'oncogene', type = 'not')
    )

HD_fisher_table <- rbind(LA_counts, HD_counts)

HD_fisher <- HD_fisher_table%>%
  select(not_oncogene_counts, oncogene_counts) %>%
  fisher.test()
print(paste("HD:",HD_fisher$p.value))