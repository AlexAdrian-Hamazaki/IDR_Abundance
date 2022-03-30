# Purpose:

# See if Highly Disordered, Low Abundant proteins are enriched for delatorious overexpression phenotypes
# compared to the proteome

# looking at Disopred3% score

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

source(file = "scripts/functions.R")


#~~~~~~~~~ Open Master ~~~~~~~~~~~~

master <- read.delim('data/merged_tables/yeast/actual_master.tsv',
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

#~~~~~~~~~ vaboori replication ~~~~~~~~~~~~
proteome <- master %>%
  filter(!is.na(disopred3perc) & !is.na(abundance))
proteome_counts <- get_int_counts(proteome, 'OE')

HD <- master  %>%
  filter(!is.na(disopred3perc) & !is.na(abundance)) %>%
  get_highest_percent(disopred3perc, disorder_percent)
HD_counts <- get_int_counts(HD, 'OE')

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

#~~~~~~~~~ Get all low abundant proteins in proteome ~~~~~~~~~~~~

LA <- master %>%
  filter(!is.na(abundance) & !is.na(disopred3perc)) %>%
  get_highest_percent(abundance, 0.33, from = 'low')

#~~~~~~~~~ Of get High Disorder low abundant proteins ~~~~~~~~~~~~

HD_LA <- master %>%
  filter(!is.na(disopred3perc) & !is.na(abundance)) %>%
  get_highest_percent(disopred3perc, disorder_percent) %>%
  get_highest_percent(abundance, abundance_percent, from = 'low')

#~~~~~~~~~ Count number of oncogenes in each  ~~~~~~~~~~~~

LA_counts <- LA %>%
  get_int_counts("OE")


HD_LA_counts <- HD_LA %>%
  get_int_counts("OE")

#~~~~~~~~~ Bring contingency table together and do fisher exact test ~~~~~~~~~~~~

fisher_table <- rbind(LA_counts, HD_LA_counts)
fisher_table <- select(fisher_table, c(counts, not_counts))
rownames(fisher_table) <- c('LA', 'HD_LA')

fisher <- fisher.test(fisher_table)
print(paste("P_value for LA_HD vs LA from the bottom = ",abundance_percent, "disorder percent", disorder_percent, "=", fisher$p.value))

#~~~~~~~~~ Now as a control see if the high abundance high disopred is enriched with OE ~~~~~~~~~~

HD_HA <- master %>%
  filter(!is.na(disopred3perc) & !is.na(abundance)) %>%
  get_highest_percent(disopred3perc, disorder_percent) %>%
  get_highest_percent(abundance, abundance_percent)

#~~~~~~~~~ Count number of oncogenes ~~~~~~~~~~~~

HD_HA_counts <- HD_HA %>%
    summarize(df_rows = nrow(HD_HA),
            OE_counts = sum(OE, na.rm = TRUE),
            not_OE_counts = sum(is.na(OE))
    )

#~~~~~~~~~ Bring contingency table together and do fisher exact test ~~~~~~~~~~~~

fisher_table_HA <- rbind(LA_counts, HD_HA_counts)
fisher_table_HA <- select(fisher_table_HA, c(OE_counts, not_OE_counts))
rownames(fisher_table_HA) <- c('LA', 'HD_HA')

fisher <- fisher.test(fisher_table_HA)
print(paste("P_value for HA_HD vs LA percent FROM THE BOTTOM = ",abundance_percent, "disorder percent", disorder_percent, "=", fisher$p.value))



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
            OE_counts = sum(OE, na.rm = TRUE),
            not_OE_counts = sum(is.na(OE))
    )

HD_HA_counts <- HD_HA %>%
    summarize(df_rows = nrow(HD_HA),
            OE_counts = sum(OE, na.rm = TRUE),
            not_OE_counts = sum(is.na(OE))
    )

# Bind contingency table

total_fisher_table <- rbind(LA_counts, HD_LA_counts, HD_MA_counts, HD_HA_counts)
rownames(total_fisher_table) <- c('LA', 'HD_LA', 'HD_MA', 'HD_HA')

#  Fishers for High abundance, mid abundance and low abundance

HD_HA_fisher <- total_fisher_table %>%
  select(OE_counts, not_OE_counts) %>%
  filter(rownames(total_fisher_table)%in% c('LA', 'HD_HA'))%>%
  fisher.test()
print(paste("HD_HA:", HD_HA_fisher$p.value))

HD_MA_fisher <- total_fisher_table %>%
  select(OE_counts, not_OE_counts) %>%
  filter(rownames(total_fisher_table)%in% c('LA', 'HD_MA'))%>%
  fisher.test()
print(paste("HD_MA:",HD_MA_fisher$p.value))


HD_LA_fisher <- total_fisher_table %>%
  select(OE_counts, not_OE_counts) %>%
  filter(rownames(total_fisher_table)%in% c('LA', 'HD_LA'))%>%
  fisher.test()
print(paste("HD_LA vs LA:",HD_LA_fisher$p.value))


# Do a fisher for the entirey of High Disorder vs Proteome

HD <- master %>%
  filter(!is.na(disopred3perc) & !is.na(abundance)) %>%
  get_highest_percent(disopred3perc, disorder_percent)

HD_counts <- HD%>%
    summarize(df_rows = nrow(HD),
            OE_counts = sum(OE, na.rm = TRUE),
            not_OE_counts = sum(is.na(OE))
    )

proteome <- master %>%
  filter(!is.na(disopred3perc) & !is.na(abundance))
proteome_counts <- proteome %>%
    summarize(df_rows = nrow(proteome),
            OE_counts = sum(OE, na.rm = TRUE),
            not_OE_counts = sum(is.na(OE))
    )

HD_fisher_table <- rbind(proteome_counts, HD_counts)

HD_fisher <- HD_fisher_table%>%
  select(OE_counts, not_OE_counts) %>%
  fisher.test()
print(paste("HD vs Proteome:",HD_fisher$p.value))