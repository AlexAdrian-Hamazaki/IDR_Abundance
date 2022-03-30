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
master$dr_all_phase_sep <- unlist(lapply(master$dr_all_phase_sep, NA_to_0))

master$dr_all_phase_sep <- factor(master$dr_all_phase_sep, levels = c('0', '1'))
levels(master$dr_all_phase_sep)


ggplot(master) +
  geom_histogram(mapping = aes(x = abundance, fill = dr_all_phase_sep), bin = 30) +
  geom_vline(mapping = aes(xintercept = median(abundance))) +
  geom_vline(mapping = aes(xintercept = quantile(master$abundance)[[4]]))

### violin plot

master <- master %>%
  group_by(dr_all_phase_sep)


ggplot(master, mapping = aes( x = dr_all_phase_sep, y = abundance) ) +
  geom_violin(fill = 'steelblue')+
  geom_boxplot(width = 0.03)+
  labs(title = "Abundance Distribution of Proteome vs Scaffolds") +
  ylab("Abundance (log10)")+
  xlab('') +
  scale_x_discrete(labels = c('Proteome', 'Scaffolds'))+
  theme(axis.text.x = element_text(face = "plain", size = 12))


#### 2 Sample T Test

non_client <- master %>%
  filter(dr_all_phase_sep == 0)

clients <- master %>%
  filter(dr_all_phase_sep == 1)

t_test <- t.test(non_client$abundance, clients$abundance)
t_test$p.value
t_test
#is a t-test ok? - no, the log10 abundance is not normal

wilcox <- wilcox.test(non_client$abundance, clients$abundance)
wilcox
