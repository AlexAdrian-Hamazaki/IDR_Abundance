# Function

# For a given Master data frame, determine the enrichment of HI or OE
# Do so by binning abundance and disorderd proteins, and calculating enrichments.


library(tidyr)
library(dplyr)
library(pheatmap)
library('RColorBrewer')
library(stringr)

source(file = 'scripts/functions.R')


get_bins <- function(abundance_groups, disorder_groups) {
  bins <- list()
  for ( i in 1:abundance_groups) {
    for ( j in 1:disorder_groups) {
      bins <- append(bins, list(paste0(i,",",j)))
    }
  }
  return(bins)
}
create_list <- function(abundance_groups, disorder_groups) {
  #Create a list for the binning of genes.
  #X axis of the matrix will be binsf rom abundance values
  #Y axis of matrix will be bins from disorder values
  #X and Y should be numbers

  bins <- get_bins(abundance_groups, disorder_groups)
  vector <- rep(NA, abundance_groups*disorder_groups)
  names(vector) <- bins
  return (vector)
}
iterate_bounds <- function(group_size, unit) {
  #get group bounds as a vector
  output_vector <- (1)
  for (i in 1:(group_size)) {
    next_bound <- ceiling(i*unit)
    output_vector <- append(output_vector, next_bound)
  }
  return (output_vector)
}
get_bounds <- function(data_frame, group_col, group_size){
  #Order your group col by descending values
  #Evenly devide your group by the the group size
  #Return the bounds of those groups as a list?

  #group_col is the column name of your target group as a string

  group <- na.omit(data_frame[ , group_col])
  len <- length(group)
  unit <- len/(group_size)

  bounds <- iterate_bounds(group_size, unit)

  return(bounds)
}
get_binned_dfs <- function(df, bounds, bounds_column){
  #Returns data frames split into different leves based on bounds, and the bounds column name
  #bounds_column is a string of the column name that you are splitting
  master_ordered <- arrange(df, !!sym(bounds_column))
  master_splits <- list()
  for (i in 1:(length(bounds)-1)) {
    if (i == 1) {
      my_master <- master_ordered[ bounds[i] : bounds[i+1], ]
      master_splits <- append(master_splits, list(my_master))
  } else {
    correct_index <- bounds[i] + 1
    my_master <- master_ordered[correct_index:bounds[i+1], ]
    master_splits <- append(master_splits, list(my_master))
    }
  }
  return(master_splits)
}
add_bin_tag <- function(master_splits, newcolname) {
  #add a new column and tag every df with "1" in that column
  #newcolname should be the column name you want to add as a str
  for ( i in 1:length(master_splits) ) {
    master_splits[[i]] <- mutate(master_splits[[i]], "{newcolname}" := i)
  }
  return(master_splits)
}
bind_master_splits <- function(master_splits) {
  # rbind your master splits
  master_new <- bind_rows(master_splits)
  return(master_new)
}
add_bin_info <- function(df, bounds, bounds_column, newcolname) {
  # df = master data frame
  # Bounds = the bounds you will split your frame at
  # bounds_column = the column name (string) that you will be using to split via your bounds
  # newcolname = the new column name you will add to indicate which bin a row is in (str)

  # This function takes in your master data frame, arranges it in descending order via a column, splits the df via
  # bounds, and then onto each of the data frames it adds a new column indicating which bin/bounds it is in
  # Then it merges the data frames back together and returns the frame.

  splits <- get_binned_dfs(df = df, bounds = bounds, bounds_column = bounds_column)
  abff <- add_bin_tag(splits, newcolname = newcolname)
  master_new <- bind_master_splits(abff)
  return(master_new)
}
place_rows_in_bins <- function(both_binned, my_list) {
  for (row_number in 1:nrow(both_binned)) {
    message(paste('Evaluating row number:', row_number))
    row <- both_binned[row_number,]
    binval <- paste0(row$abundance_bin,",",row$disorder_bin)
    for ( listname in names(my_list)) {
        if ( binval == listname & !is.na(my_list[binval]) ) {
          my_list[[binval]]<- rbind(my_list[[binval]], row)
          break
        } else if ( binval == listname & is.na(my_list[binval]) ) {
          my_list[binval]<- list(row)
          break
        }
    }

  }
  return (my_list)
}
get_p_value <- function(enrichment_df, alternative = 'greater') {
  # For use in an lapply
  # get a p-value for each df in enrichment_list2
  fisher_table <- enrichment_df %>%
    select(in_bin, not_in_bin)
  fisher_test <- fisher.test(fisher_table, alternative = alternative)
  p_val <- fisher_test$p.value
  return (p_val)
}
get_effect_size <- function(enrichment_df) {
  # For use in an lapply
  # get a effect size for each df in enrichment_list1:
  # Effect size will be the percentage of genes that are enriched for whatever you are looking at (HI or OE typically)
  percentage <- enrichment_df[[1,1]] / enrichment_df [[1,2]]

  return(percentage*100)
}
make_matrix <- function(abundance_groups, disorder_groups) {
  #Make a matrix for which you will put P-values in
  my_matrix <- matrix(NA, ncol = abundance_groups, nrow = disorder_groups)
  return(my_matrix)
}
fill_matrix <- function(list, my_matrix) {
  #Take the list of our P-values and make it into a X by Y matrix
  for (name in names(list)) {
    name_split <- str_split(name, ',')
    name_int <- lapply(name_split, as.integer)
    abundance_bin <- name_int[[1]][1]
    disorder_bin <- name_int[[1]][2]
    my_matrix[disorder_bin, abundance_bin] <- list[[name]]
  }
  rownames(my_matrix) <- c(paste0("disorder",1:nrow(my_matrix)))
  colnames(my_matrix) <- c(paste0("abundance",1:ncol(my_matrix)))
  return(my_matrix)
}
make_breaks <- function(p_values, alpha = 0.05, gradients = 5) {
  #make breaks for pheatmap
  dist <- unlist(p_values[p_values<=0.05])
  probabilities <- seq(from = 0, to = 1, length.out = gradients)
  quantile <- quantile(dist, probs = probabilities)
  quantile <- append(quantile, alpha)
  quantile <- append(quantile, 1)
  return(quantile)
}
make_colors <- function(gradients = 4){
  #make colors for pheatmap
  colors <- rev(brewer.pal( n = (gradients), "Reds"))
  colors <- append(colors, 'grey')
  return(colors)
}
create_bins <- function(master, abundance_groups, disorder_groups) {
      master <- master %>%
      filter(!is.na(disopred3perc) & !is.na(abundance))

    abundance_bounds <- get_bounds(data_frame = master, group_col = 'abundance', group_size = abundance_groups)
    disorder_bounds <- get_bounds(data_frame = master, group_col = 'disopred3perc', group_size = disorder_groups)
    abundance_binned <- add_bin_info(df = master,
                                   bounds = abundance_bounds,
                                   bounds_column = 'abundance',
                                   newcolname = 'abundance_bin')
    both_binned <- add_bin_info(df=abundance_binned,
                              bounds = disorder_bounds,
                              bounds_column = 'disopred3perc',
                              newcolname = 'disorder_bin')

    my_list <- create_list(abundance_groups, disorder_groups)
    binned_genes <- place_rows_in_bins(both_binned, my_list)
  return(binned_genes)
}
make_and_fill_matrix <- function(values_list, abundance_groups, disorder_groups) {
  #Make and fill your p_value matrix for the purpose of putting it in a pheatmap
  # Values list can be a list of p-values or effect_sizes values
  value_matrix <- make_matrix(abundance_groups = abundance_groups,
                            disorder_groups = disorder_groups)
  value_matrix <- fill_matrix(values_list, value_matrix)
  value_matrix <- value_matrix[nrow(value_matrix):1 , ]
  return(value_matrix)
}

################## Functions specific to this script
get_not_hits <- function(enrichment_data_frame, total_hits) {
  # for use in an lapply loop
  # for each data frame in list, append the number of "not hits"
  # Return data frame with appended column of not hits
  hits <- enrichment_data_frame[1, 'counts']
  not_hits <- total_hits - hits
  enrichment_data_frame[1,'not_in_bin'] <-  not_hits
  colnames(enrichment_data_frame) <- c('total', 'in_bin', 'not_hits','not_in_bin')


  return(enrichment_data_frame)
}
append_not_hits <- function(enrichment_list) {
  enrichment_table <- do.call(rbind, enrichment_list)
  total_hits <- sum(enrichment_table$counts)
  enrichment_list2 <- lapply(enrichment_list, get_not_hits, total_hits)
  return(enrichment_list2)
}
get_proteome_bin_total <- function(enrichment_list) {
  # for use in an lapply
  # extracts the value for total
  return (as.data.frame(enrichment_list[1 , 'total']))
}
append_not_in_bin <- function(proteome_bin_table, total_proteins) {
  # for use in lapply loop
  # Append the total number of proteins not in each bin
  not_in_bin <- total_proteins - proteome_bin_table[1,1]
  proteome_bin_table <- proteome_bin_table %>%
    mutate(not_in_bin = not_in_bin)
  colnames(proteome_bin_table) <- c('in_bin', 'not_in_bin')
  return(proteome_bin_table)
}
append_proteome_to_enrichment <- function(enrichment_list, not_hits) {
  # append the nth table in not_hits to the first table in enrichment list
  fisher_tables <- list()
  names <- names(enrichment_list)
  for (index in 1:length(enrichment_list)) {
    table <- rbind(enrichment_list[[index]], not_hits[[index]])
    rownames(table) <- c('hit', 'proteome')

    fisher_tables <- append(fisher_tables, list(table))
  }
  names(fisher_tables) <- names
  return(fisher_tables)
}
##################3 Values to edit
stopifnot(FALSE)

abundance_groups <-  5
disorder_groups <-  5
gradients <-  5
alpha <-  0.01
enrichment_column <- 'is_HI'
disorder_column <- 'aa'
abundance_column <- 'ff'

master <- read.delim('data/merged_tables/human/actual_master.tsv',
                     sep = "\t")
# log2+1 abundance
master$abundance <- log2(master$abundance+1)


binned_genes <- create_bins(master, abundance_groups, disorder_groups)
enrichment_list <- lapply(binned_genes, get_int_counts, sym(enrichment_column))
enrichment_list <- append_not_hits(enrichment_list)

proteome_bin_total <- lapply(enrichment_list, get_proteome_bin_total)
total_proteins <- do.call(sum, proteome_bin_total)
not_hits <- lapply(proteome_bin_total, append_not_in_bin, total_proteins)

enrichment_list <- lapply(enrichment_list, select, c(2,4))


fisher_tables <- append_proteome_to_enrichment(enrichment_list, not_hits)

p_values_list <- lapply(fisher_tables, get_p_value)
effect_sizes_list <- lapply(fisher_tables, get_effect_size)

# p_breaks <- make_breaks(unlist(p_value_matrix), alpha = alpha, gradients = gradients)
# p_colors <- make_colors(gradients = gradients)
p_value_matrix <- make_and_fill_matrix(p_values_list, abundance_groups, disorder_groups)
pheatmap(mat = p_value_matrix,
         main = paste0("Abundance:",abundance_groups, " x ", "Disorder:",disorder_groups, " Pvalue heatmap for GROUP enrichment in HI Proteins"),
          breaks = c(0, alpha, 1),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = c('red4', 'grey')
)

effect_size_matrix <- make_and_fill_matrix(effect_sizes_list, abundance_groups, disorder_groups)

pheatmap(mat = effect_size_matrix,
         main = paste0("Abundance:",abundance_groups, " x ", "Disorder:",disorder_groups, " % Enriched Heatmap for GROUP enrichment in HI"),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = brewer.pal(n = 9, name = 'Reds'),
         legend_labels = c("%")
)
