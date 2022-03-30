#Purpose:

#  Do Fisher Exact Tests with baseline as proteome high abundance
#  Compare Highly disordered disopred vs baseline


library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)

source('scripts/functions.R')

master <- read.delim('data/merged_tables/yeast/actual_master.tsv',
                     sep = "\t")
# log2+1 abundance
master$abundance <- log2(master$abundance+1)


#~~~~~~~~~ Global Parameters ~~~~~~~~~~~~
abundance_percent <- 0.4
disorder_percent <- 0.333333
# For human_total, this top percent will be taken as "highly abundant" and "highly disordered"
print(paste("Abundance Percent:", abundance_percent))
print(paste("Disorder Percent:", disorder_percent))


#~~~~~~~~~ Get all high abundant proteins in proteome ~~~~~~~~~~~~
HA <- master %>%
  filter(!is.na(abundance) & !is.na(disopred3perc)) %>%
  get_highest_percent(abundance, abundance_percent)

HA_counts <- get_int_counts(HA, column_name =  "HI")

#~~~~~~~~~ Get  High abundant, High disorder proteins~~~~~~~~~~~~
HD_HA <- master %>%
  filter(!is.na(abundance) & !is.na(disopred3perc)) %>%
  get_highest_percent(disopred3perc, disorder_percent) %>%
  get_highest_percent(abundance, abundance_percent)

HD_HA_counts <- get_int_counts(HD_HA, column_name = "HI")

#~~~~~~~~~ Get counts for the mid abundant, High disorder proteins~~~~~~~~~~~~
HD_MA <- master %>%
  filter(!is.na(abundance) & !is.na(disopred3perc)) %>%
  get_highest_percent(disopred3perc, disorder_percent)%>%
  get_middle_percent(abundance, abundance_percent)

HD_MA_counts <-  get_int_counts(HD_MA, column_name = "HI")

#~~~~~~~~~ Get counts for the low abundant, High disorder proteins~~~~~~~~~~~~
HD_LA <- master %>%
    filter(!is.na(abundance) & !is.na(disopred3perc)) %>%
      get_highest_percent(disopred3perc, disorder_percent)%>%
      get_lowest_percent(abundance, abundance_percent)

HD_LA_counts <- get_int_counts(HD_LA, column_name = "HI")

#~~~~~~~~~ bind together contingency table ~~~~~~~~~~~~

fisher_total <- rbind(HA_counts, HD_HA_counts, HD_MA_counts, HD_LA_counts)
rownames(fisher_total) <- c('HA', 'HD_HA', 'HD_MA', 'HD_LA')

#~~~~~~~~~ Do Fisher tests for X vs HA~~~~~~~~~~~~
HD_HA_fisher <- fisher_total %>%
  filter(rownames(fisher_total) %in% c('HA', 'HD_HA'))%>%
  fisher.test()
print(paste("P.val HD_HA vs HA:",HD_HA_fisher$p.value))

HD_MA_fisher <- fisher_total %>%
  filter(rownames(fisher_total) %in% c('HA', 'HD_MA'))%>%
  fisher.test()
print(paste("P.val HD_MA vs HA:",HD_MA_fisher$p.value))

HD_LA_fisher <- fisher_total %>%
  filter(rownames(fisher_total) %in% c('HA', 'HD_LA'))%>%
  fisher.test()
print(paste("P.val HD_LA vs HA:" , HD_LA_fisher$p.value))

#~~~~~~~~~ Graphing ~~~~~~~~~~~~

fisher_elongate <- fisher_total %>%
  select(-c(total)) %>%
  mutate(names  = rownames(fisher_total)) %>%
    melt()
name_levels <- c('HD_HA', 'HD_MA', 'HD_LA', "HA")
ggplot(data = fisher_elongate) +
  geom_bar(mapping =  aes(x = factor(names, levels = name_levels), y = value, fill = variable),
           position='stack', stat = 'identity')+
  labs(title = "HI counts for HD proteins across different abundance levels",
       x = "Different Abundance Levels",
       y = "counts",
      fill = "HI Protein Presence")+
  scale_fill_discrete(labels = c('HI Proteins', 'Not HI Proteins'))+
  scale_x_discrete(labels = name_levels)

stopifnot(FALSE)


#~~~~~~~~~ Get counts for the High Abundant, Mid disorder proteins~~~~~~~~~~~~
MD_HA <- master %>%
    filter(!is.na(abundance) & !is.na(disopred3perc)) %>%
      get_middle_percent(disopred3perc, disorder_percent)%>%
      get_highest_percent(abundance, abundance_percent)
#~~~~~~~~~ Get counts for HI genes in Low abundant, High disorder proteins~~~~~~~~~~~~

MD_HA_counts <- get_int_counts(MD_HA, column_name = "HI")

#~~~~~~~~~ Get counts for the High Abundant, Low disorder proteins~~~~~~~~~~~~
LD_HA <- master %>%
    filter(!is.na(abundance) & !is.na(disopred3perc)) %>%
      get_lowest_percent(disopred3perc, disorder_percent)%>%
      get_highest_percent(abundance, abundance_percent)
#~~~~~~~~~ Get counts for HI genes in Low abundant, High disorder proteins~~~~~~~~~~~~

LD_HA_counts <- get_int_counts(LD_HA, column_name = "HI")


#~~~~~~~~~ Get all high disorder proteins in proteome ~~~~~~~~~~~~#

HD <- master %>%
  filter(!is.na(abundance) & !is.na(disopred3perc)) %>%
  get_highest_percent(disopred3perc, disorder_percent)

HD_counts <- get_int_counts(HD, column_name = "HI")
#HD counts should contain it the total amount of HI genes we can get in the HD_HA+HD_MA+HD_LA
#~~~~~~~~~ bind together contingency table ~~~~~~~~~~~~

fisher_total <- rbind(HA_counts,
                      HD_counts,
                      HD_HA_counts,
                      HD_MA_counts,
                      HD_LA_counts,
                      MD_HA_counts,
                      LD_HA_counts
)

rownames(fisher_total) <-  as.character(expression(HA_counts,
                      HD_counts,
                      HD_HA_counts,
                      HD_MA_counts,
                      HD_LA_counts,
                      MD_HA_counts,
                      LD_HA_counts))

#~~~~~~~~~ Do Fisher tests ~~~~~~~~~~~~
HD_HA_fisher <- fisher_total %>%
  filter(rownames(fisher_total) %in% c('HA_counts', 'HD_HA_counts')) %>%
  fisher.test()
print(paste("HD_HA vs HA:",HD_HA_fisher$p.value))

HD_MA_fisher <- fisher_total %>%
  filter(rownames(fisher_total) %in% c('HA_counts', 'HD_MA_counts'))%>%
  fisher.test()
print(paste("HD_MA vs HA:",HD_MA_fisher$p.value))

HD_LA_fisher <- fisher_total %>%
  filter(rownames(fisher_total) %in% c('HA_counts', 'HD_LA_counts'))%>%
  fisher.test()
print(paste("HD_LA vs HA:" , HD_LA_fisher$p.value))

MD_HA_fisher <- fisher_total %>%
  filter(rownames(fisher_total) %in% c('HD_HA_counts', 'MD_HA_counts')) %>%
  fisher.test()
print(paste("HD_HA vs MD_HA:" , HD_LA_fisher$p.value))
#Seems actually like there's a de-enrichment of HI at HD_HA

stopifnot(FALSE)

#~~~~~~~~~ Do a fisher test for the entirey of High Disorder ~~~~~~~~~~~~
HD <- master %>%
    filter(!is.na(abundance) & !is.na(disopred3perc)) %>%
    get_highest_percent(disopred3perc, disorder_percent)
HD_count <- summarize_HI(HD)

HD_fisher_table <- rbind(HA_counts, HD_count)
HD_fisher <- fisher.test(HD_fisher_table)
print(paste("P_value for TOTAL HD at abundance percentage of  = ",abundance_percent, "disorder percent", disorder_percent, "=", HD_fisher$p.value))

#~~~~~~~~~ Do a fisher test for the High disorder, LOW abundant but from the bottom % ~~~~~~~~~~~~
#This will select our lowest abundance bin from the bottom instead of getting it by splitting the mid-low bins
#This is just to ensure that any p-values from splitting the mid-low bins are not due to small sample size of those bins
HD_LA_BOT <- master %>%
    filter(!is.na(abundance) & !is.na(disopred3perc)) %>%
    get_highest_percent(disopred3perc, disorder_percent) %>%
    get_highest_percent(abundance, abundance_percent, from = 'low')
HD_LA_BOT_counts <- summarize_HI(HD_LA_BOT)
HD_LA_BOT_fisher_table <- rbind(HA_counts, HD_LA_BOT_counts)
HD_LA_BOT_fisher <- fisher.test(HD_LA_BOT_fisher_table)
print(paste("P_value for HD_HA FROM THE BOTTOM at abundance percentage of  = ",abundance_percent, "disorder percent", disorder_percent, "=", HD_LA_BOT_fisher$p.value))

#~~~~~~~~~ Global Parameters ~~~~~~~~~~~~
abundance_percent <- 0.4
disorder_percent <- 0.333333
# For human_total, this top percent will be taken as "highly abundant" and "highly disordered"
print(paste("Abundance Percent:", abundance_percent))
print(paste("Disorder Percent:", disorder_percent))



stopifnot(FALSE)


####

##ARCHIVE

####



human_master <- read.delim(file = "data/merged_tables/human/human_master.tsv",
                           sep = "\t")
#Get High disorder, high abundance genes, and then count the number of HI's
HD_MA <- human_master %>%
  get_highest_percent(disopred, disorder_percent) %>%
  get_highest_percent(abundance, abundance_percent) %>%
  summarize(HI_counts = count_string(HI_bin, "HI"),
            nHI_counts = count_string(HI_bin, "nHI"),
            miss_counts = count_string(HI_bin, "missing"),
            combined_misses = nHI_counts + miss_counts)
rownames(HD_HA) <- "HD_HA"

# Get high abundant genes, and then count the number of HI's
HA <- human_master %>%
  get_highest_percent(abundance, abundance_percent) %>%
  summarize(HI_counts = count_HI(HI_bin, ""),
            nHI_counts = count_HI(HI_bin, "nHI"),
            miss_counts = count_HI(HI_bin, "missing"),
            combined_misses = nHI_counts + miss_counts)
rownames(HA) <- "HA"


# Build the fisher exact test table
contingency_table<- rbind(HA, HD_HA)


#~~~~~~~~~ In medium disordered: Fisher Exact Test using entire proteome as comparison~~~~~~~~~~~~


MD_HA <- human_master %>%
  get_middle_percent(disopred, disorder_percent) %>%
  get_highest_percent(abundance, abundance_percent) %>%
  summarize(HI_counts = count_HI(HI_bin, ""),
          nHI_counts = count_HI(HI_bin, "nHI"),
          miss_counts = count_HI(HI_bin, "missing"),
          combined_misses = nHI_counts + miss_counts)
rownames(MD_HA) <- "MD_HA"

#we already have human _high abundance data as HA

# Build the fisher exact test table
contingency_table<- rbind(contingency_table, MD_HA)


#~~~~~~~~~ In low disordered: Fisher Exact Test using entire proteome as comparison~~~~~~~~~~~~

LD_HA <- human_master %>%
  get_lowest_percent(disopred, disorder_percent) %>%
  get_highest_percent(abundance, abundance_percent) %>%
  summarize(HI_counts = count_HI(HI_bin, ""),
          nHI_counts = count_HI(HI_bin, "nHI"),
          miss_counts = count_HI(HI_bin, "missing"),
          combined_misses = nHI_counts + miss_counts)
rownames(LD_HA) <- "LD_HA"

#we already have human _high abundance data as HA

contingency_table<- rbind(contingency_table, LD_HA)

#~~~~~~~~~ Doing Fisher Stats ~~~~~~~~~~~~

# HD_HA vs HA Fisher
HD_fisher_table <- contingency_table %>%
  filter(rownames(contingency_table) %in% c('HD_HA', 'HA')) %>%
  select(HI_counts, combined_misses)

HD_fisher_test <-  fisher.test(HD_fisher_table)
print(paste("HD_fisher p-value : ", HD_fisher_test$p.value))

# MD_HA vs HA Fisher
MD_fisher_table <- contingency_table %>%
  filter(rownames(contingency_table) %in% c('MD_HA', 'HA')) %>%
  select(HI_counts, combined_misses)

MD_fisher_test <-  fisher.test(MD_fisher_table)
print(paste("MD_fisher p-value : ", MD_fisher_test$p.value))

# LD_HA vs HA Fisher
LD_fisher_table <- contingency_table %>%
  filter(rownames(contingency_table) %in% c('LD_HA', 'HA')) %>%
  select(HI_counts, combined_misses)

LD_fisher_test <-  fisher.test(LD_fisher_table)
print(paste("LD_fisher p-value : ", LD_fisher_test$p.value))

stopifnot(FALSE)
#~~~~~~~~~ Graphing ~~~~~~~~~~~~

contingency_table <-  contingency_table %>%
  mutate("Percent" = HI_counts / combined_misses * 100)
percent_diffs_gg <- ggplot(contingency_table) +
  geom_col(mapping = aes( x = factor(rownames(contingency_table),
                                     levels = rownames(contingency_table)),
                          y = Percent),
           fill = 'steelblue') +
  theme_classic()+
  labs(title = "Percent of Happloinsifficient Genes",
       x= "Disorder levels",
       y = "Percent Happloinsuficiency") +
  scale_x_discrete(labels = c("HighAbunProteome", "HighAbun", "MidAbun", "LowAbun"))
percent_diffs_gg

ggsave(percent_diffs_gg, file = "../../../../../figures/Happloinsufficiency/HI_at_diff_abun_levels/abun_0.4.png")