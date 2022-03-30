#Purpose:
#Identify if proteins that are SG are enriched for high disorder high abundance
#In humans
#Here, the background we are using is the proteome

library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)


source('scripts/binning_functions.R')
source('scripts/enrichment_summary_functions.R')

master <- read.delim('data/merged_tables/human/actual_master.tsv',
                     sep = "\t")
nsg <- read.delim(file = 'data/phase_sep/Erich_nsg/nsg_human.txt', header = FALSE)
#~~~~~~~~~ Global Parameters ~~~~~~~~~~~~

disorder_col <-  'disopred3perc'
abundance_col <-  'abundance'
enrichment_column <-  'stress_granule'
abundance_groups <-  3
disorder_groups <-  3
alpha = 0.01

### 1) Get Stress Granule Proteins

stress_granule_proteins <- filter(master, uniprot_id %in% nsg$V1) %>%
  mutate(stress_granule = as.integer(1))
not_stress_granule_proteins <- filter(master, !uniprot_id %in% nsg$V1) %>%
  mutate(stress_granule = as.integer(0))

### Append to master if the protein is a stress granule or not
# if it is, then the 'stress_granule' column will be == 1

master <- rbind(stress_granule_proteins, not_stress_granule_proteins)



# log2+1 abundance
master$abundance <- log2(master$abundance+1)


### remove NAs

master <- filter(master,
                 !is.na(!!sym(abundance_col)) & !is.na(!!sym(disorder_col)))


### 2) Bin out Genes for Disorder and Abundance

# Create bins
binned_genes <- create_bins(master,
                            abundance_groups = abundance_groups,
                            disorder_groups = disorder_groups,
                            disorder_col = disorder_col,
                            abundance_col = abundance_col)

######## Fisher Exact Test to determine enrichment

### Getting proteome summary data:
# For each bin, how many proteins are in that bin, and how many are not in that bin


proteome_diso_abun_list <- lapply(binned_genes, get_proteome_in_bin, master)

### 3) Get total amount of SG genes in each bin, and how many are in total


enrichment_in_bins <- lapply(binned_genes, get_enrichment_in_bin, proteome = master, enrichment_column_name = enrichment_column)

### 4) Add the proteome diso_abun_lists onto each of the bins in enrichment_list


enrichment_list <- rbind_dfs_in_list(enrichment_in_bins, proteome_diso_abun_list)


### 5) Get Effect Size for each of the bins

#Appnds effect size straight onto enrichment_list

enrichment_list <- lapply(enrichment_list, get_effect_size, hit_column ='num_in_bin', miss_column = 'num_not_in_bin')

### 6) Get P-values for each bin

enrichment_p_values <- lapply(enrichment_list, get_p_value, alpha = alpha, hit_column ='num_in_bin', miss_column = 'num_not_in_bin')

### 7) Visualize Effect Size

# To vizuaize, you need to transfer your values into a heatmap for pheatmap

effect_sizes <- unlist(lapply(enrichment_list, pull_effect_sizes))

effect_size_matrix <- matrix(effect_sizes*100,
                             nrow =  disorder_groups,
                             ncol = abundance_groups)

# Flips the matrix over its horizontal plane so that the highest disorder highest abundance bin is top right
effect_size_matrix <- ((effect_size_matrix [ c(nrow(effect_size_matrix) : 1) , ]))


### Visualize the effect size
pheatmap(mat = effect_size_matrix,
         main = paste0("Abundance:",abundance_groups, " x ", "Disorder:",disorder_groups, " % Enriched Heatmap"),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = brewer.pal(n = 9, name = 'Reds'),
         legend = TRUE,
         annotation_legend = TRUE,
         legend_labels = c(5,1,3,3,6),
         show_rownames = TRUE,
         labels_col = 'Abundance Bins',
         labels_row = 'Disorder Bins',
         angle_col = 0,
         angle_row = 90
)


### 8) Visualize P-Values



p_values_double <- unlist(enrichment_p_values)
p_value_matrix <- matrix(p_values_double, ncol = abundance_groups, nrow = disorder_groups)

# Flips the matrix over its horizontal plane so that the highest disorder highest abundance bin is top right
p_value_matrix <- ((p_value_matrix[ c(nrow(p_value_matrix) : 1) , ]))

### Visualize

pheatmap(mat = p_value_matrix,
         main = paste0("Abundance:",abundance_groups, " x ", "Disorder:",disorder_groups, " Pvalue heatmap for enrichment of stress granules"),
         breaks = c(0, alpha, 1),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = c('red4', 'grey'),
         labels_col = 'Abundance Bins',
         labels_row = 'Disorder Bins',
         angle_col = 0,
         angle_row = 90
)
