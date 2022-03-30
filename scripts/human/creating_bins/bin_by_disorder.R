
# Takes in human_total_data.tsv Arranges via disopred% score, and then adds a disopred_bin column, seperating into low, mid and high
# Then does the same for iupred% score

# Then adds a dose_bin: dose scores of 1,2,3 indicate HI, 40 or 30 indicate not HI, and 0 Indicates no data.
# Then re-saves human_total_data.tsv with the new data. But only saves uniprot_id and the bin.
# The intention of the save is to just join it to your human_total later on


library(tidyr)
library(dplyr)
library(ggplot2)

human_total <- read.delim('data/merged_tables/human/human_total_data.tsv')

# Take log of abundance
human_total$abundance <- log2(human_total$abundance+1)

#~~~~~~~~~ SELECT TOP BIN % ~~~~~~~~~~~~~
# If you want to increase/decrease your top bin size you  can just change this parameter to change it them
top_bin_percent <- 0.35

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



#~~~~~~~~~ Create Disopred% bins ~~~~~~~~~~~~~
disopred_bins <- human_total %>%
  filter(!is.na(disopred)) %>%
  arrange(disopred)%>%
  select(uniprot_id, disopred)


bin_index <- creat_bins(disopred_bins, top_bin_percent)

disopred_bins <- disopred_bins %>%
  mutate(
         disopred_percent_bin = c(rep('low', bin_index[1]),
                                  rep('mid', bin_index[2]-bin_index[1]),
                                  rep('high', bin_index[3]-bin_index[2])))

disopred_low <- disopred_bins[1:bin_index[1], ]
disopred_mid <- disopred_bins[(bin_index[1]+1):bin_index[2], ]
disopred_high <- disopred_bins[(bin_index[2]+1):bin_index[3], ]

ggplot(disopred_low) +
  geom_histogram(mapping = aes(disopred)) +
  labs(title = "disopred%score_low")
ggplot(disopred_mid) +
  geom_histogram(mapping = aes(disopred)) +
  labs(title = "disopred%score_mid")
ggplot(disopred_high) +
  geom_histogram(mapping = aes(disopred)) +
  labs(title = "disopred%score_high")

disopred_bins <- select(disopred_bins, -disopred)
human_total <- left_join(human_total, disopred_bins, by = "uniprot_id", .keepall = FALSE)
stopifnot(sum(is.na(human_total$disopred_percent_bin)) ==sum(is.na(human_total$disopred)))

#~~~~~~~~~ Create iupred% bins ~~~~~~~~~~~~~
iupred_bins <- human_total %>%
  filter(!is.na(iupred)) %>%
  arrange(iupred)%>%
  select(uniprot_id, iupred)

bin_index <- creat_bins(iupred_bins, top_bin_percent)

iupred_bins <- iupred_bins %>%
  mutate(
         iupred_percent_bin = c(rep('low', bin_index[1]),
                                  rep('mid', bin_index[2]-bin_index[1]),
                                  rep('high', bin_index[3]-bin_index[2])))

iupred_low <- iupred_bins[1:bin_index[1], ]
iupred_mid <- iupred_bins[(bin_index[1]+1):bin_index[2], ]
iupred_high <- iupred_bins[(bin_index[2]+1):bin_index[3], ]

ggplot(iupred_low) +
  geom_histogram(mapping = aes(iupred)) +
  labs(title = "iupred%score_low")
ggplot(iupred_mid) +
  geom_histogram(mapping = aes(iupred)) +
  labs(title = "iupred%score_mid")
ggplot(iupred_high) +
  geom_histogram(mapping = aes(iupred)) +
  labs(title = "iupred%score_high")

iupred_bins <- select(iupred_bins, -iupred)
human_total <- left_join(human_total, iupred_bins, by = "uniprot_id", .keepall = FALSE)
stopifnot(sum(is.na(human_total$iupred_percent_bin)) ==sum(is.na(human_total$iupred)))

human_total <- human_total %>%
  select(uniprot_id, disopred_percent_bin, iupred_percent_bin)



# Write table contining the disopred percent bins, iupred percent bins, and a uniprot ID so we can merge
# with human_total later
dir.create('data/merged_tables/human', showWarnings = FALSE)
write.table(x = human_total, file = "data/merged_tables/human/human_disorder_bins.tsv", sep = '\t')
