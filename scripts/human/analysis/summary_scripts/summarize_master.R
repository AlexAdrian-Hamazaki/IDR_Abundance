library(tidyverse)

master <- read.delim(file = 'data/merged_tables/human/actual_master.tsv', sep = "\t")



get_filled <- function(df, colname){
  data <- df %>%
    filter(!is.na(!!sym(colname)))

  return (nrow(data))
}
get_filled(master, 'is_HI')


summary <- data.frame('Abundance' = get_filled(master, 'abundance'),
                      'Disorder' = get_filled(master, 'disopred3perc'),
                      'Haploinsufficient' = get_filled(master, 'is_HI'),
                      'Oncogene' = get_filled(master, 'oncogene'),
                      'PS Driver'= get_filled(master, 'ps_high_conf'),
                      'PS Client' = get_filled(master, 'drpls_clients'))

master_yeast <- read.delim(file = 'data/merged_tables/yeast/actual_master.tsv', sep = "\t")

summary_yeast <- data.frame('Abundance' = get_filled(master_yeast, 'abundance'),
                      'Disorder' = get_filled(master_yeast, 'disopred3perc'),

                      'Haploinsufficient' = get_filled(master_yeast, 'HI'),
                      'OPsapko' = get_filled(master_yeast, 'OE'),

                      'OPmorril' = get_filled(master_yeast, 'new_oe'),

                      'PS Driver'= get_filled(master_yeast, 'drlps_drivers'),
                      'PS Client' = get_filled(master_yeast, 'drpls_clients'))


master_human <- master %>%
  filter(!is.na(abundance))%>%
  filter(!is.na(disopred3perc))

master_yeast_2 <- master_yeast %>%
  filter(!is.na(abundance))%>%
  filter(!is.na(disopred3perc))



