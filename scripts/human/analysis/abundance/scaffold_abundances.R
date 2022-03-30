library(tidyverse)

master <- read.delim(file = "data/merged_tables/human/actual_master.tsv", sep = "\t")

master$abundance <- log10(master$abundance+1)

master <- master %>%
  filter(!is.na(abundance))

NA_to_0 <- function(optional) {
  if (is.na(optional)) {
    return (0)
  } else {
    return (1)
  }
}
master$ps_high_conf <- unlist(lapply(master$ps_high_conf, NA_to_0))

master$ps_high_conf <- factor(master$ps_high_conf, levels = c('0', '1'))
levels(master$ps_high_conf)


ggplot(master) +
  geom_histogram(mapping = aes(x = abundance, fill = ps_high_conf), bin = 30) +
  geom_vline(mapping = aes(xintercept = median(abundance))) +
  geom_vline(mapping = aes(xintercept = quantile(master$abundance)[[4]]))

### violin plot

master <- master %>%
  group_by(ps_high_conf)


ggplot(master, mapping = aes( x = ps_high_conf, y = abundance) ) +
  geom_violin(fill = 'steelblue')+
  geom_boxplot(width = 0.03)+
  labs(title = "Abundance Distribution of Proteome vs Scaffolds") +
  ylab("Abundance (log10)")+
  xlab('') +
  scale_x_discrete(labels = c('Proteome', 'Scaffolds'))+
  theme(axis.text.x = element_text(face = "plain", size = 12))


#### 2 Sample T Test

non_scaf <- master %>%
  filter(ps_high_conf == 0)

scaf <- master %>%
  filter(ps_high_conf == 1)

t_test <- t.test(non_scaf$abundance, scaf$abundance)
t_test$p.value

#is a t-test ok? - no, the log10 abundance is not normal

wilcox <- wilcox.test(non_scaf$abundance, scaf$abundance)
wilcox$p.value