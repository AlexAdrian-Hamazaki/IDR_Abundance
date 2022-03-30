library(tidyverse)
library(reshape2)
female <- data.frame(val = c( 1.9, 1.6, 1.4, 1.1, 1.6, 1.8, 1.9, 1.7, 1.5, 1.8, 1.7, 1.7,
1.8, 1.7, 1.8, 2.0, 1.8, 1.7, 1.6, 1.6, 1.5), cat = 'fem')

intact_male <- data.frame(val = c( 1.9, 1.2, 1.0, 0.9, 1.4, 1.0, 1.3, 1.4, 1.1, 1.0, 1.4,
1.2, 1.4, 1.4, 1.5, 1.5, 1.1, 1.4, 1.3, 1.3, 1.3), cat = 'in_m')

removed_minor_male <- data.frame( val = c(1.2, 1.0, 0.9, 0.8, 1.2, 0.9, 1.1, 1.1, 1.3, 1.3,
1.3, 1.1, 1.4, 1.5, 1.4, 1.4, 1.2, 1.4, 1.3, 1.2, 1.4), cat = 'rem_min_m')

removed_major_male <- data.frame(val = c(1.2, 0.9, 1.4, 1.2, 1.2, 1.6, 1.9, 1.4, 1.4, 1.4,
1.6, 1.4, 1.7, 1.3, 1.5, 1.2, 1.3, 1.6, 1.5, 1.5, 1.5), cat = 'rem_maj_m')

crabs <- rbind(female, intact_male, removed_minor_male, removed_major_male)


my_lm <- lm(val~cat, data = crabs)

my_anova <- anova(my_lm)
my_anova