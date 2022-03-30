# Take the total human data that has bins for low,mid and high disorder.
# For each group, bin by abundance scores into low, mid and high.
# Save a data frame for the low, mid, and high disorder.

# Do this for diopred3%, iupred%, and later disorder metrics as well

library(dplyr)
library(tidyr)
library(ggplot2)

#~~~~~~~~~ Load master table ~~~~~~~~~~~~~
human_total <- read.delim(file = 'data/merged_tables/human/human_total_data.tsv', sep ='\t')
#~~~~~~~~~ Add on disorder information ~~~~~~~~~~~~~
human_disorder_bins <- read.delim(file = "data/merged_tables/human/human_disorder_bins.tsv", sep = "\t")
human_total <- left_join(human_total, human_disorder_bins, by = "uniprot_id")
#~~~~~~~~~ Add on dose sensitive information ~~~~~~~~~~~~~
human_HI_bins <- read.delim(file = "data/merged_tables/human/human_HI_bins.tsv", sep = "\t")
human_total <- left_join(human_total, human_HI_bins, by = "uniprot_id")

#Just take log2 abundance right away
human_total$abundance <- log2(human_total$abundance+1)

#~~~~~~~~~ SELECT TOP BIN % ~~~~~~~~~~~~~
# If you want to increase/decrease your top bin size you  can just change this parameter to change it them
top_bin_percent <- 0.65

creat_bins <- function(dataframe, top_bin_percent){
  # Given a data frame, split the data frame into low, mid and high
  # High will be designated first, and it will be the size of the top_bin_percent parameter
  # Mid and low will be split in half from the remaining rows
  # Return a tuple of length 3, Indicating the slices for the table
  # @PARAM: top_bin_percent: must be in 0.33 format. not 33%
  rows <- nrow(dataframe)
  top_percent <- floor(rows*top_bin_percent+0.5)
  second_index <- rows - top_percent
  first_index <- floor(second_index/2+0.5)

  return(c(first_index, second_index, rows))
}

#~~~~~~~~~~~~~~~DISOPRED3PERCENT~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~DISOPRED3PERCENT~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~DISOPRED3PERCENT~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~ Create high disopred% abundance bins ~~~~~~~~~~~~~
disopred_per_high <- human_total %>%
  filter(disopred_percent_bin == "high")%>%
  filter(!is.na(abundance)) %>%
  arrange(abundance)


bin_index <- creat_bins(disopred_per_high, top_bin_percent)

disopred_per_high <-
  mutate(disopred_per_high, abundance_bin =
    c(rep('low', bin_index[1]),
      rep('mid', bin_index[2]-bin_index[1]),
      rep('high', bin_index[3]-bin_index[2])))

DH_abun_low <- disopred_per_high[1:bin_index[1], ]
DH_abun_mid <- disopred_per_high[(bin_index[1]+1):bin_index[2], ]
DH_abun_high <- disopred_per_high[(bin_index[2]+1):bin_index[3], ]

ggplot(DH_abun_low) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of High Disopred%, low Abundance", x = "log2abun")
ggplot(DH_abun_mid) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of High Disopred%, mid Abundance",x = "log2abun")
ggplot(DH_abun_high) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of High Disopred%, high Abundance", x = "log2abun")

#~~~~~~~~~ Create mid disopred% abundance bins ~~~~~~~~~~~~~
disopred_per_mid <- human_total %>%
  filter(disopred_percent_bin == "mid")%>%
  filter(!is.na(abundance)) %>%
  arrange(abundance)

bin_index <- creat_bins(disopred_per_mid, top_bin_percent)


disopred_per_mid <-
  mutate(disopred_per_mid, abundance_bin =
    c(rep('low', bin_index[1]),
      rep('mid', bin_index[2]-bin_index[1]),
      rep('high', bin_index[3]-bin_index[2])))

DM_abun_low <- disopred_per_mid[1:bin_index[1], ]
DM_abun_mid <- disopred_per_mid[(bin_index[1]+1):bin_index[2], ]
DM_abun_high <- disopred_per_mid[(bin_index[2]+1):bin_index[3], ]

ggplot(DM_abun_low) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of mid Disopred%, low Abundance", x = "log2abun")
ggplot(DM_abun_mid) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of mid Disopred%, mid Abundance",x = "log2abun")
ggplot(DM_abun_high) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of mid Disopred%, high Abundance", x = "log2abun")

#~~~~~~~~~ Create low disopred% abundance bins ~~~~~~~~~~~~~
disopred_per_low <- human_total %>%
  filter(disopred_percent_bin == "low")%>%
  filter(!is.na(abundance)) %>%
  arrange(abundance)

bin_index <- creat_bins(disopred_per_low, top_bin_percent)


disopred_per_low <-
  mutate(disopred_per_low, abundance_bin =
    c(rep('low', bin_index[1]),
      rep('mid', bin_index[2]-bin_index[1]),
      rep('high', bin_index[3]-bin_index[2])))

DL_abun_low <- disopred_per_low[1:bin_index[1], ]
DL_abun_mid <- disopred_per_low[(bin_index[1]+1):bin_index[2], ]
DL_abun_high <- disopred_per_low[(bin_index[2]+1):bin_index[3], ]

ggplot(DL_abun_low) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of low Disopred%, low Abundance", x = "log2abun")
ggplot(DL_abun_mid) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of low Disopred%, mid Abundance",x = "log2abun")
ggplot(DL_abun_high) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of low Disopred%, high Abundance", x = "log2abun")

#~~~~~~~~~ save data frames  ~~~~~~~~~~~~~
# 1 data frame is saved for low disopred, mid, and high disopred.
# Each data frame has abundance binned
dir.create("data/merged_tables/human/disopred3%bins", showWarnings = FALSE)
write.table(x = disopred_per_high,
            file = "data/merged_tables/human/disopred3%bins/high_disopred.tsv",
            sep = "\t")
write.table(x = disopred_per_mid,
            file = "data/merged_tables/human/disopred3%bins/mid_disopred.tsv",
            sep = "\t")
write.table(x = disopred_per_low,
            file = "data/merged_tables/human/disopred3%bins/low_disopred.tsv",
            sep = "\t")

#~~~~~~~~~~~~~~~IUPREDPERCENT~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~IUPREDPERCENT~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~IUPREDPERCENT~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~ Create high iupred% abundance bins ~~~~~~~~~~~~~
iupred_per_high <- human_total %>%
  filter(iupred_percent_bin == "high")%>%
  filter(!is.na(abundance)) %>%
  arrange(abundance)

bin_index <- creat_bins(iupred_per_high, top_bin_percent)


iupred_per_high <-
  mutate(iupred_per_high, abundance_bin =
    c(rep('low', bin_index[1]),
      rep('mid', bin_index[2]-bin_index[1]),
      rep('high', bin_index[3]-bin_index[2])))

IH_abun_low <- iupred_per_high[1:bin_index[1], ]
IH_abun_mid <- iupred_per_high[(bin_index[1]+1):bin_index[2], ]
IH_abun_high <- iupred_per_high[(bin_index[2]+1):bin_index[3], ]

ggplot(IH_abun_low) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of High Iupred%, low Abundance", x = "log2abun")
ggplot(IH_abun_mid) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of High Iupred%, mid Abundance",x = "log2abun")
ggplot(IH_abun_high) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of High Iupred%, high Abundance", x = "log2abun")

#~~~~~~~~~ Create mid iupred% abundance bins ~~~~~~~~~~~~~
iupred_per_mid <- human_total %>%
  filter(iupred_percent_bin == "mid")%>%
  filter(!is.na(abundance)) %>%
  arrange(abundance)

bin_index <- creat_bins(iupred_per_mid, top_bin_percent)


iupred_per_mid <-
  mutate(iupred_per_mid, abundance_bin =
    c(rep('low', bin_index[1]),
      rep('mid', bin_index[2]-bin_index[1]),
      rep('high', bin_index[3]-bin_index[2])))

IM_abun_low <- iupred_per_mid[1:bin_index[1], ]
IM_abun_mid <- iupred_per_mid[(bin_index[1]+1):bin_index[2], ]
IM_abun_high <- iupred_per_mid[(bin_index[2]+1):bin_index[3], ]

ggplot(IM_abun_low) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of mid Iupred%, low Abundance", x = "log2abun")
ggplot(IM_abun_mid) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of mid Iupred%, mid Abundance",x = "log2abun")
ggplot(IM_abun_high) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of mid Iupred%, high Abundance", x = "log2abun")

#~~~~~~~~~ Create low iupred% abundance bins ~~~~~~~~~~~~~
iupred_per_low <- human_total %>%
  filter(iupred_percent_bin == "low")%>%
  filter(!is.na(abundance)) %>%
  arrange(abundance)

bin_index <- creat_bins(iupred_per_low, top_bin_percent)


iupred_per_low <-
  mutate(iupred_per_low, abundance_bin =
    c(rep('low', bin_index[1]),
      rep('mid', bin_index[2]-bin_index[1]),
      rep('high', bin_index[3]-bin_index[2])))

IL_abun_low <- iupred_per_low[1:bin_index[1], ]
IL_abun_mid <- iupred_per_low[(bin_index[1]+1):bin_index[2], ]
IL_abun_high <- iupred_per_low[(bin_index[2]+1):bin_index[3], ]

ggplot(IL_abun_low) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of low Iupred%, low Abundance", x = "log2abun")
ggplot(IL_abun_mid) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of low Iupred%, mid Abundance",x = "log2abun")
ggplot(IL_abun_high) +
  geom_histogram( mapping = aes(abundance))+
  labs(title="Dispersion of low Iupred%, high Abundance", x = "log2abun")

#~~~~~~~~~ save data frames  ~~~~~~~~~~~~~
# 1 data frame is saved for low iupred, mid, and high iupred.
# Each data frame has abundance binned
dir.create("data/merged_tables/human/iupred%bins", showWarnings = FALSE)
write.table(x = iupred_per_high,
            file = "data/merged_tables/human/iupred%bins/high_iupred.tsv",
            sep = "\t")
write.table(x = iupred_per_mid,
            file = "data/merged_tables/human/iupred%bins/mid_iupred.tsv",
            sep = "\t")
write.table(x = iupred_per_low,
            file = "data/merged_tables/human/iupred%bins/low_iupred.tsv",
            sep = "\t")

message("Ran successfully")
stopifnot(FALSE)






# HI data for highly abundant
dose_disopred_high_abun_high <- disopred_high_abun_high %>%
  filter(!is.na(HI_desc))
dose_hits_high <- dose_disopred_high_abun_high %>%
  filter(HI_score == 1 | HI_score ==2 |HI_score == 3 )
dose_misses <- dose_disopred_high_abun_high %>%
  filter(HI_score == 30 | HI_score ==40)
dose_nodata <- dose_disopred_high_abun_high %>%
  filter(HI_score == 0)
nrow(dose_hits_high) + nrow(dose_misses) + nrow(dose_nodata) == nrow(dose_disopred_high_abun_high)

print(paste(nrow(dose_hits_high) , 'were identified as high disopred, high abundance, and HI'))
print(paste(nrow(dose_misses) , 'were identified as high disopred, high abundance, and low HI'))
print(paste(nrow(dose_nodata) , 'were identified as high disopred, high abundance, and  no data for HI'))
print(nrow(dose_hits_high)/nrow(dose_misses))

write.table(x = dose_hits_high, file = 'data/merged_tables/split_by_disopred/high_disopred_splits/high_abundance/hits')
write.table(x = dose_misses, file = 'data/merged_tables/split_by_disopred/high_disopred_splits/high_abundance/misses')
write.table(x = dose_nodata, file = 'data/merged_tables/split_by_disopred/high_disopred_splits/high_abundance/nodata')


# HI data for lowly abundant
dose_disopred_high_abun_low <- disopred_high_abun_low %>%
  filter(!is.na(HI_desc))
dose_hits_low <- dose_disopred_high_abun_low %>%
  filter(HI_score == 1 | HI_score ==2 |HI_score == 3 )
dose_misses <- dose_disopred_high_abun_low %>%
  filter(HI_score == 30 | HI_score ==40)
dose_nodata <- dose_disopred_high_abun_low %>%
  filter(HI_score == 0)
nrow(dose_hits) + nrow(dose_misses) + nrow(dose_nodata) == nrow(dose_disopred_high_abun_low)

write.table(x = dose_hits_low, file = 'data/merged_tables/split_by_disopred/high_disopred_splits/low_abundance/hits')
write.table(x = dose_misses, file = 'data/merged_tables/split_by_disopred/high_disopred_splits/low_abundance/misses')
write.table(x = dose_nodata, file = 'data/merged_tables/split_by_disopred/high_disopred_splits/low_abundance/nodata')

print(paste(nrow(dose_hits) , 'were identified as high disopred, high abundance, and HI'))
print(paste(nrow(dose_misses) , 'were identified as high disopred, high abundance, and low HI'))
print(paste(nrow(dose_nodata) , 'were identified as high disopred, high abundance, and  no data for HI'))
print(nrow(dose_hits)/nrow(dose_misses))
# HI data for mid abundant
dose_disopred_high_abun_mid <- disopred_high_abun_mid %>%
  filter(!is.na(HI_desc))
dose_hits_mid <- dose_disopred_high_abun_mid %>%
  filter(HI_score == 1 | HI_score ==2 |HI_score == 3 )
dose_misses <- dose_disopred_high_abun_mid %>%
  filter(HI_score == 30 | HI_score ==40)
dose_nodata <- dose_disopred_high_abun_mid %>%
  filter(HI_score == 0)
nrow(dose_hits) + nrow(dose_misses) + nrow(dose_nodata) == nrow(dose_disopred_high_abun_mid)

write.table(x = dose_hits_mid, file = 'data/merged_tables/split_by_disopred/high_disopred_splits/mid_abundance/hits')
write.table(x = dose_misses, file = 'data/merged_tables/split_by_disopred/high_disopred_splits/mid_abundance/misses')
write.table(x = dose_nodata, file = 'data/merged_tables/split_by_disopred/high_disopred_splits/mid_abundance/nodata')

print(paste(nrow(dose_hits_mid) , 'were identified as high disopred, high abundance, and HI'))
print(paste(nrow(dose_misses) , 'were identified as high disopred, high abundance, and low HI'))
print(paste(nrow(dose_nodata) , 'were identified as high disopred, high abundance, and  no data for HI'))
print(nrow(dose_hits)/nrow(dose_misses))



#TODO bootstrap null?

