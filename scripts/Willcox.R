
library(tidyr)
library(dplyr)
library(ggplot2)

source(file = 'scripts/functions.R')

master <- read.delim('data/merged_tables/yeast/actual_master.tsv',
                     sep = "\t")
master <- master %>%
  filter(!is.na(abundance) & !is.na(disopred3perc))

master$abundance <- log2(master$abundance+1)



target_column <- 'HI'
# For Humans use 'is_HI' for yeast use 'HI'

abundance_column <- 'abundance'
disorder_column <- 'disopred3perc'

master <- mutate_at(master, target_column, ~replace(., is.na(.), 0))
master$HI <- as.factor(as.integer(master$HI))

enrichments <- master %>%
  group_by(!!sym(target_column)) %>%
  group_split()

hits <- enrichments[[2]]
not_hits <- enrichments[[1]]

ggplot(hits)+
  geom_histogram(mapping = aes(abundance), fill = 'steelblue', binwidth = 0.5) +
  labs(title =(paste('Abundance distribution for', target_column)),
       x= 'Log2(Abundance+1)')
ggplot(not_hits)+
  geom_histogram(mapping = aes(abundance), fill = 'steelblue', binwidth = 0.5) +
  labs(title =(paste('Abundance distribution for not', target_column)),
       x = 'Log2(Abundance+1)')

ggplot(hits) +
  geom_histogram(mapping = aes(!!sym(disorder_column)),
                 fill = 'steelblue',
                 binwidth = 5)+
  labs( title = paste(disorder_column, 'Disorder distribution for', target_column))

ggplot(not_hits) +
  geom_histogram(mapping = aes(!!sym(disorder_column)),
                 fill = 'steelblue',
                 binwidth = 5)+
  labs( title = paste(disorder_column, 'Disorder distribution for not', target_column))

ggplot(master, mapping = aes(x = !!sym(target_column),  y = abundance)) +
  geom_boxplot() +
  labs(title = paste("Abundance differences for", target_column), "and not") +
  xlab ('Hit and not Hit groups') +
  ylab( "Log2(Abundance+1)") +
  scale_x_discrete(labels = c('Not Hits', target_column))

ggplot(master, mapping = aes(x = !!sym(target_column),  y = disopred3perc)) +
  geom_boxplot() +
  labs(title = paste("Disorder differences for", target_column, "and not")) +
  xlab ('Hit and not Hit groups') +
  scale_x_discrete(labels = c('Not Hits', target_column))

abundance_wilcox <- wilcox.test(hits$abundance, not_hits$abundance, alternative = 'greater')
disorder_wilcox <- wilcox.test(hits$disopred3perc, not_hits$disopred3perc)

message(paste('Abundance Wilcox P value:', abundance_wilcox$p.value))
message(paste('Disorder Wilcox P value:', disorder_wilcox$p.value))

