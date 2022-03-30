#Purpose:
#Identify if proteins that are disordered and abundant are enriched for haploinsufficiency
#In human
#Here, the background we are using is the proteome

library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


source('scripts/binning_functions.R')
source('scripts/enrichment_summary_functions.R')

master <- read.delim('data/merged_tables/human/actual_master.tsv',
                     sep = "\t")
#~~~~~~~~~ Global Parameters ~~~~~~~~~~~~

disorder_col <-  'disopred3perc'
abundance_col <-  'abundance'
enrichment_column <-  'is_HI'
abundance_groups <-  5
disorder_groups <-  5
alpha = 0.01




# log2+1 abundance
master$abundance <- log2(master$abundance+1)


### 1) remove NAs

master <- filter(master,
                 !is.na(!!sym(abundance_col)) & !is.na(!!sym(disorder_col)))


### 2) Bin out Genes for Disorder and Abundance
# Create bins
binned_genes <- create_bins(master,
                            abundance_groups = abundance_groups,
                            disorder_groups = disorder_groups,
                            disorder_col = disorder_col,
                            abundance_col = abundance_col)

### Fisher Exact Test to determine if Stress Granule proteins are enriched

enrichment_list <- lapply(binned_genes, summarize_enrichment, enrichment_column)

### 3) Get the summary enrichment for the proteome as our comparison

proteome_enrichment <- summarize_enrichment(master, enrichment_column)

### 4) Add the proteome enrichment summary table onto each of the bins in enrichment_list

enrichment_list <- lapply(enrichment_list, rbind, proteome_enrichment)

### 5) Get Effect Size for each of the bins

#Appnds effect size straight onto enrichment_list
enrichment_list <- lapply(enrichment_list, get_effect_size)


### 6) Get P-values for each bin

enrichment_p_values <- lapply(enrichment_list, get_p_value, alpha = alpha)

### 7) Visualize Effect Size

# To vizuaize, you need to transfer your values into a heatmap for pheatmap

effect_sizes <- unlist(lapply(enrichment_list, pull_effect_sizes))

effect_size_matrix <- matrix(effect_sizes*100,
                             nrow =  disorder_groups,
                             ncol = abundance_groups)

# Flips the matrix over its horizontal plane so that the highest disorder highest abundance bin is top right
effect_size_matrix <- ((effect_size_matrix [ c(nrow(effect_size_matrix) : 1) , ]))


### Visualize the effect size
stopifnot(FALSE)
#png(filename = 'figures/human/HI/5x5HI.png', width = 2000, height = 950, res = 300)
pheatmap(mat = effect_size_matrix,
         main = paste0("Human: Percent Enrichment of Haploinsufficiency"),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = brewer.pal(n = 9, name = 'Reds'),
         legend = TRUE,
         annotation_legend = TRUE,
         legend_labels = c(5,1,3,3,6),
         show_rownames = TRUE,
         labels_row = '                 ',
         angle_col = 0,
         angle_row = 90
)
#dev.off()


### 8) Visualize P-Values



p_values_double <- unlist(enrichment_p_values)
p_value_matrix <- matrix(p_values_double, ncol = abundance_groups, nrow = disorder_groups)

# Flips the matrix over its horizontal plane so that the highest disorder highest abundance bin is top right
p_value_matrix <- ((p_value_matrix[ c(nrow(p_value_matrix) : 1) , ]))

### Visualize
png(filename = 'figures/human/HI/5x5HI_p.png', width = 2000, height = 950, res = 300)



pheatmap(mat = p_value_matrix,
         main = paste0("Human: Percent Enrichment of Haploinsufficiency P Values"),
         breaks = c(0, alpha, 1),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = c('red4', 'grey'),
         labels_col = 'Abundance Bins',
         labels_row = 'Disorder Bins',
         angle_col = 0,
         angle_row = 90
)
dev.off()