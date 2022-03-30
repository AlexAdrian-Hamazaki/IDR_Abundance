library(tidyverse)

master <- read.delim(file = 'data/merged_tables/yeast/actual_master.tsv')

#Filter for HI and PS
ps_col <- 'dr_all_phase_sep'

HI <- master %>%
  select(HI, sym(ps_col)) %>%
  filter(HI == 1)

sum(HI$dr_all_phase_sep, na.rm = TRUE)
