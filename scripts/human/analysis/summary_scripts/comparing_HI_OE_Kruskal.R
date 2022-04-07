library(tidyverse)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(googlesheets4)


master <- read.delim('data/merged_tables/human/actual_master.tsv',
                     sep = "\t")
#~~~~~~~~~ Global Parameters ~~~~~~~~~~~~





oncogenes <- master %>%
  filter(oncogene == 1)

HI <- master %>%
  filter(is_HI == 1)

proteome <- master %>%
  filter(!uniprot_id %in% oncogenes$uniprot_id) %>%
  filter(!uniprot_id %in% HI$uniprot_id)

proteome[,'bin'] <- 'Proteome'
oncogenes[,'bin'] <- 'Oncogenes'
HI[,'bin'] <- 'Haploinsufficient'

merged <- rbind(proteome, oncogenes)
merged <- rbind(merged, HI)
merged$bin <- factor(merged$bin, levels = c('Proteome', 'Oncogenes', 'Haploinsufficient'))


ggplot(merged, mapping = aes(x = bin, y = disopred3perc) ) +
  geom_violin(fill = 'steelblue')+
  geom_boxplot(width = 0.03)+
  labs(title =
         "Yeast: Disorder Distributions of
     Phase Separating Proteins
          ") +
  ylab("% Disorder")+
  theme_classic()+
  xlab(' ') +
  theme(axis.text.x = element_text(face = "plain", size = 12), plot.title = element_text(face = "bold"))
ggsave(filename = 'figures/yeast/comparisons/yeast_disorder.png', dpi = 300, width = 1200, height = 1300, units = 'px')


###

master <- read.delim(file ='data/merged_tables/yeast/actual_master.tsv', sep = "\t")

HI <- master %>%
  filter(HI == 1)
sum(HI$OE, na.rm= TRUE)


?kruskal.test()

