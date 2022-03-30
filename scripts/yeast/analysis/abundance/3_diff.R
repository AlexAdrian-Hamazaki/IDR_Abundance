library(tidyverse)
library(dunn.test)

master <- read.delim(file = "data/merged_tables/yeast/actual_master.tsv", sep = "\t")


# Process Abundance
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

# Factorize the client phase separating proteins
master$drpls_clients <- unlist(lapply(master$drpls_clients, NA_to_0))

master$drpls_clients <- factor(master$drpls_clients, levels = c('0', '1'))
levels(master$drpls_clients)


# Factorize the driver proteins
master$drlps_drivers <- unlist(lapply(master$drlps_drivers, NA_to_0))

master$drlps_drivers <- factor(master$drlps_drivers, levels = c('0', '1'))
levels(master$drlps_drivers)


# Get scaffolds
scafs <- master %>%
  filter(drlps_drivers == 1)

# Get clients
clients <- master %>%
  filter(drpls_clients == 1)

# Get proteome background
proteome <- master %>%
  filter(drpls_clients != 1)

#Get how many scafs are also somehow clients
stranges <- clients %>%
  filter(clients$uniprot_id %in% scafs$uniprot_id)

#Remove the client+drivers from the client frame. They are high confidence drivers
clients <- clients %>%
  filter(! clients$uniprot_id %in% stranges$uniprot_id)


# Add factors for merging of dfs for plotting later
proteome$bin <- factor('Proteome')
clients$bin <- factor('Client')

scafs$bin <- factor('Driver')


#merge 3 tables
merged <- rbind(proteome, clients)
merged <- rbind(merged, scafs)


#anova (kruskal wallis test)
krus <- kruskal.test(abundance ~ bin, merged)
krus

#dunn

dunn <- dunn.test(x = merged$abundance, g = merged$bin, method = 'bh')
dunn$P.adjusted
dunn

# Violin Plot
stopifnot(FALSE)
ggplot(merged, mapping = aes(x = bin, y = abundance) ) +
  geom_violin(fill = 'steelblue')+
  geom_boxplot(width = 0.03)+
  labs(title = "Human: Abundance Distributions of Phase Separating Proteins") +
  ylab("Abundance (log10)")+
  scale_x_discrete(labels = c('Proteome', 'Clients', 'Drivers'))+
  scale_y_continuous(limits = c(0, 5)) +
  xlab(' ') +
  theme(axis.text.x = element_text(face = "plain", size = 12), plot.title = element_text(face = "bold"))


ggsave(filename = "figures/yeast/abundances/abundances_PS.png", width = 6.5, height = 3, dpi = 300)


