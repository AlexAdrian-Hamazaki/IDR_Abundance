#Purpose:
#Identify if proteins that are disordered and abundant are enriched for dose sensitivity
#In Yeast
#Here, the background we are using is the proteome

library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


source('scripts/binning_functions.R')
source('scripts/enrichment_summary_functions.R')


master <- readRDS(file = 'data/merged_tables/yeast/actual_master.rds')
#~~~~~~~~~ Global Parameters ~~~~~~~~~~~~

disorder_col <-  'disopred3perc'
abundance_col <-  'abundance'
enrichment_column <-  'OE'
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

######## Fisher Exact Test to determine enrichment ########

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
png(filename = 'figures/yeast/OP/5x5newOE.png', width = 2000, height = 850, res = 300)

a <- pheatmap(mat = effect_size_matrix,
         main = paste0("Yeast: Percent Enrichment of Overexpression Phenotypes"),
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
print(a)
dev.off()




### 8) Visualize P-Values



p_values_double <- unlist(enrichment_p_values)
p_value_matrix <- matrix(p_values_double, ncol = abundance_groups, nrow = disorder_groups)

# Flips the matrix over its horizontal plane so that the highest disorder highest abundance bin is top right
p_value_matrix <- ((p_value_matrix[ c(nrow(p_value_matrix) : 1) , ]))

### Visualize
stopifnot(FALSE)

png(filename = 'figures/yeast/OP/5x5newOE_p.png', width = 2000, height = 850, res = 300)

# pheatmap(mat = p_value_matrix,
#          main = paste0("Yeast: Percent Enrichment of Overexpression Phenotypes P Values"),
#          breaks = c(log2(0), log2(alpha), log2(1)),
#          cluster_cols = FALSE,
#          cluster_rows = FALSE,
#          color = c('red4', 'grey'),
#          labels_row = '                 ',
#          angle_col = 0,
#          angle_row = 90
# )


pheatmap(mat = p_value_matrix,
         main = paste0("Yeast: Percent Enrichment of Overexpression Phenotypes P Values"),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = rev(brewer.pal(n = 9, name = 'Reds')),
         labels_col = 'Abundance Bins',
         labels_row = 'Disorder Bins',
         angle_col = 0,
         angle_row = 90
)

dev.off()

png(filename = 'figures/yeast/OP/5x5oldOE_p.png', width = 1400, height = 850, res = 300)

color <- brewer.pal(n = 5, name = 'Reds')
color <- rev(append(color, 'black'))
pheatmap(mat = p_value_matrix,
         main = paste0("Yeast: Sapko OP Enrichment P Values"),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = color,
         labels_row = '         ',
         angle_col = 0,
         angle_row = 90,
         breaks = c(0, alpha, 0.1, 0.25, 0.5, 0.75, 1),
)





dev.off()

ff <- format(p_value_matrix, scientific = TRUE, trim = TRUE, digits = 3)


write.csv(ff, file = 'figures/test/p_value_matrix.csv')