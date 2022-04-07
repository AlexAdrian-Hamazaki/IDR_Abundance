library(tidyverse)
library(ggplot2)
library(dunn.test)


### YEAST DISORDER
master <- read.delim(file = 'data/merged_tables/yeast/actual_master.tsv', sep = "\t")

drivers <- master %>%
  filter(drlps_drivers == 1)

clients <- master %>%
  filter(drpls_clients == 1)

proteome <- master %>%
  filter(!uniprot_id %in% drivers$uniprot_id) %>%
  filter(!uniprot_id %in% clients$uniprot_id)

proteome[,'bin'] <- 'Proteome'
drivers[,'bin'] <- 'Drivers'
clients[,'bin'] <- 'Clients'

merged <- rbind(proteome, clients)
merged <- rbind(merged, drivers)
merged$bin <- factor(merged$bin, levels = c('Proteome', 'Clients', 'Drivers'))


ggplot(merged, mapping = aes(x = bin, y = disopred3perc) ) +
  geom_violin(fill = 'steelblue')+
  geom_boxplot(width = 0.03)+
  labs(title =
         "Yeast: Disorder Distributions of
     Phase Separating Proteins
          ") +
  ylab("% Disorder")+
  theme_classic()+
  scale_x_discrete(labels = c('Proteome', 'Clients', 'Drivers'))+
  xlab(' ') +
  theme(axis.text.x = element_text(face = "plain", size = 12), plot.title = element_text(face = "bold"))
ggsave(filename = 'figures/yeast/comparisons/yeast_disorder.png', dpi = 300, width = 1200, height = 1300, units = 'px')

#anova (kruskal wallis test)
krus <- kruskal.test(disopred3perc ~ bin, merged)
krus

#dunn

dunn <- dunn.test(x = merged$disopred3perc, g = merged$bin, method = 'bh')
dunn


### HUMAN DISORDER #!!!!!!!~~~~~~~~~~~~~~~~~~~~
master <- read.delim(file = 'data/merged_tables/human/actual_master.tsv', sep = "\t")

drivers <- master %>%
  filter(ps_high_conf== 1)

clients <- master %>%
  filter(drpls_clients == 1)

proteome <- master %>%
  filter(!uniprot_id %in% drivers$uniprot_id) %>%
  filter(!uniprot_id %in% clients$uniprot_id)

proteome[,'bin'] <- 'Proteome'
drivers[,'bin'] <- 'Drivers'
clients[,'bin'] <- 'Clients'

merged <- rbind(proteome, clients)
merged <- rbind(merged, drivers)
merged$bin <- factor(merged$bin, levels = c('Proteome', 'Clients', 'Drivers'))


ggplot(merged, mapping = aes(x = bin, y = disopred3perc) ) +
  geom_violin(fill = 'steelblue')+
  geom_boxplot(width = 0.03)+
  labs(title =
         "Human: Disorder Distributions of
       Phase Separating Proteins
          ") +
  ylab("% Disorder")+
  theme_classic()+
  scale_x_discrete(labels = c('Proteome', 'Clients', 'Drivers'))+
  xlab(' ') +
  theme(axis.text.x = element_text(face = "plain", size = 12), plot.title = element_text(face = "bold"))
ggsave(filename = 'figures/human/comparisons/human_disorder.png', dpi = 300, width = 1200, height = 1300, units = 'px')

#anova (kruskal wallis test)
krus <- kruskal.test(disopred3perc ~ bin, merged)
krus

#dunn

dunn <- dunn.test(x = merged$disopred3perc, g = merged$bin, method = 'bh')
dunn



#####~~~~~~~ HUMAN ABUNDANCE
master <- read.delim(file = 'data/merged_tables/human/actual_master.tsv', sep = "\t")

drivers <- master %>%
  filter(ps_high_conf== 1)

clients <- master %>%
  filter(drpls_clients == 1)

proteome <- master %>%
  filter(!uniprot_id %in% drivers$uniprot_id) %>%
  filter(!uniprot_id %in% clients$uniprot_id)

proteome[,'bin'] <- 'Proteome'
drivers[,'bin'] <- 'Drivers'
clients[,'bin'] <- 'Clients'

merged <- rbind(proteome, clients)
merged <- rbind(merged, drivers)
merged$bin <- factor(merged$bin, levels = c('Proteome', 'Clients', 'Drivers'))

merged <- merged %>%
  filter(!is.na(abundance))

ggplot(merged, mapping = aes(x = bin, y = log2(abundance+1))) +
  geom_violin(fill = 'steelblue')+
  geom_boxplot(width = 0.03)+
  labs(title =
         "Human: Abundance Distributions of
       Phase Separating Proteins
          ") +
  ylab(bquote(~Log[2]~ 'Abundance (ppm)'))+
  theme_classic()+
  scale_x_discrete(labels = c('Proteome', 'Clients', 'Drivers'))+
  xlab(' ') +
  theme(axis.text.x = element_text(face = "plain", size = 12), plot.title = element_text(face = "bold"))
ggsave(filename = 'figures/human/comparisons/human_abundance.png', dpi = 300, width = 1200, height = 1300, units = 'px')

#anova (kruskal wallis test)
krus <- kruskal.test(abundance ~ bin, merged)
krus

#dunn

dunn <- dunn.test(x = merged$abundance, g = merged$bin, method = 'bh')
dunn



### YEAST ABUDANCEN~~~~~~~~~~~`
master <- read.delim(file = 'data/merged_tables/yeast/actual_master.tsv', sep = "\t")

drivers <- master %>%
  filter(drlps_drivers == 1)

clients <- master %>%
  filter(drpls_clients == 1)

proteome <- master %>%
  filter(!uniprot_id %in% drivers$uniprot_id) %>%
  filter(!uniprot_id %in% clients$uniprot_id)

proteome[,'bin'] <- 'Proteome'
drivers[,'bin'] <- 'Drivers'
clients[,'bin'] <- 'Clients'

merged <- rbind(proteome, clients)
merged <- rbind(merged, drivers)
merged$bin <- factor(merged$bin, levels = c('Proteome', 'Clients', 'Drivers'))

merged <- merged %>%
  filter(!is.na(abundance))

ggplot(merged, mapping = aes(x = bin, y = log2(abundance+1))) +
  geom_violin(fill = 'steelblue')+
  geom_boxplot(width = 0.03)+
  labs(title =
         "Yeast: Abundance Distributions of
       Phase Separating Proteins
          ") +
  ylab(bquote(~Log[2]~ 'Abundance (ppm)'))+
  theme_classic()+
  scale_x_discrete(labels = c('Proteome', 'Clients', 'Drivers'))+
  xlab(' ') +
  theme(axis.text.x = element_text(face = "plain", size = 12), plot.title = element_text(face = "bold"))
ggsave(filename = 'figures/yeast/comparisons/yeast_abundance.png', dpi = 300, width = 1200, height = 1300, units = 'px')

#anova (kruskal wallis test)
krus <- kruskal.test(abundance ~ bin, merged)
krus

#dunn

dunn <- dunn.test(x = merged$abundance, g = merged$bin, method = 'bh')
dunn