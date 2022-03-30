
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
HI_column <-  'is_HI'
phase_sep_column <- "ps_high_conf"
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

### 3) In each bin, filter for only HI proteins

filter_HI <- function(bin, column_name) {
  HI_bin <- filter(bin, !!sym(column_name) == 1)
  return(HI_bin)
}

HI_bins <- lapply(binned_genes, filter_HI, HI_column)

### 4) Calc enrichment of phase sep proteins in HI bins

enrichment_list <- lapply(HI_bins, summarize_enrichment, enrichment_column = phase_sep_column)

### 5) calc enrichment of proteome as a whole

all_HI <- master %>%
  filter(is_HI == TRUE)
proteome_enrichment <- summarize_enrichment(all_HI, enrichment_column = phase_sep_column)

### 6) Append proteome onto each bin

enrichment_list <- lapply(enrichment_list, rbind, proteome_enrichment)

### 7) Get effect size
enrichment_list <- lapply(enrichment_list, get_effect_size)

### 8) Fisher exact Test

enrichment_p_values <- lapply(enrichment_list, get_p_value, alpha = alpha)

### 9) Visualize Effect Size

# To vizuaize, you need to transfer your values into a heatmap for pheatmap


effect_sizes <- unlist(lapply(enrichment_list, pull_effect_sizes))

effect_size_matrix <- matrix(effect_sizes*100,
                             nrow =  disorder_groups,
                             ncol = abundance_groups)

# Flips the matrix over its horizontal plane so that the highest disorder highest abundance bin is top right
effect_size_matrix <- ((effect_size_matrix [ c(nrow(effect_size_matrix) : 1) , ]))

pheatmap(mat = effect_size_matrix,
         main = paste0("Human: Percent Enrichment of Phase Separation"),
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
