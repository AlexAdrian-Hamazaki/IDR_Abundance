# Function

# For a given Master data frame, determine the enrichment of HI or OE
# Do so by binning abundance and disorderd proteins, and calculating enrichments.


library(tidyr)
library(dplyr)
library(pheatmap)
library('RColorBrewer')
library(stringr)

source(file = 'scripts/functions.R')

create_matrix <- function(abundance_groups, disorder_groups) {
  #Create a matrix for the binning of genes.
  #X axis of the matrix will be binsf rom abundance values
  #Y axis of matrix will be bins from disorder values
  # X and Y should be numbers
  my_matrix <- matrix(data = NA,
                      nrow = abundance_groups,
                      ncol = disorder_groups,
                      dimnames = list(c(1:abundance_groups), c(1:disorder_groups))
  )
  return (my_matrix)
}

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
#master_splits[[i]] <- cbind(master_splits[[i]], c(rep(i, n = length(master_splits[[i]]))))

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

create_enrichment_list <- function(binned_genes){
  enrichments <- c(rep(NA, length(binned_genes)))
  names(enrichments) <- names(binned_genes)
  return(enrichments)
}

get_proteome_enrichment <- function(master, column_name) {
  # Get the enrichment of a certain value (HI or OE) for the entire proteome
  # For use of Fisher exact test later
  enrichment_analyis <- get_int_counts(master,sym(column_name))
  return (enrichment_analyis)
}

bind_proteome_enrichment <- function(df, proteome_enrichment) {
  #For use in an lapply
  # rbind the proteome_enrichment df onto each df in enrichment_list
  output <- rbind(df, proteome_enrichment)
  rownames(output) <- c('bin', 'proteome')
  return(output)
}

get_p_value <- function(enrichment_df) {
  # For use in an lapply
  # get a p-value for each df in enrichment_list2
  fisher_table <- enrichment_df %>%
    select(counts, not_counts)
  fisher_test <- fisher.test(fisher_table, alternative = 'greater')
  p_val <- fisher_test$p.value
  return (p_val)
}

get_effect_size <- function(enrichment_df) {
    # For use in an lapply
  # get a effect size for each df in enrichment_list1:
  # Effect size will be the percentage of genes that are enriched for whatever you are looking at (HI or OE typically)
  percentage <- enrichment_df[[2]] / enrichment_df [[3]]
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


# master <- master %>%
#   filter(!is.na(disopred3perc) & !is.na(abundance))
# abundance_groups <- 3
# disorder_groups <- 2
# abundance_bounds <- get_bounds(data_frame = master, group_col = 'abundance', group_size = abundance_groups)
# disorder_bounds <- get_bounds(data_frame = master, group_col = 'disopred3perc', group_size = disorder_groups)
# abundance_binned <- add_bin_info(df = master,
#                                  bounds = abundance_bounds,
#                                  bounds_column = 'abundance',
#                                  newcolname = 'abundance_bin')
# both_binned <- add_bin_info(df=abundance_binned,
#                             bounds = disorder_bounds,
#                             bounds_column = 'disopred3perc',
#                             newcolname = 'disorder_bin')
#
# my_list <- create_list(abundance_groups, disorder_groups)
# binned_genes <- place_rows_in_bins(both_binned, my_list)
# enrichment_list <- create_enrichment_list(binned_genes)
# enrichment_list <- lapply(binned_genes, get_int_counts, "HI")
# proteome_enrichment <- get_proteome_enrichment(master, 'HI')
# enrichment_list2 <- lapply(enrichment_list, bind_proteome_enrichment, proteome_enrichment)
# p_values <- lapply(enrichment_list2, get_p_value)
# p_value_matrix <- make_matrix(p_values, disorder_groups = disorder_groups, abundance_groups = abundance_groups)
# p_value_matrix <- fill_matrix(p_values, p_value_matrix)
# p_value_matrix <- p_value_matrix[nrow(p_value_matrix):1 , ]
#
# pheatmap(mat = p_value_matrix,
#          main = paste0(abundance_groups, "x", disorder_groups, " Pvalue heatmap"),
#          breaks = c(0, 0.05, 1),
#          cluster_cols = FALSE,
#          cluster_rows = FALSE,
#          color = c('red', 'grey'))


get_matrices <- function(master, abundance_groups, disorder_groups, enrichment_column) {
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
    enrichment_list <- create_enrichment_list(binned_genes)
    enrichment_list <- lapply(binned_genes, get_int_counts, sym(enrichment_column))
    proteome_enrichment <- get_proteome_enrichment(master, sym(enrichment_column))
    enrichment_list2 <- lapply(enrichment_list, bind_proteome_enrichment, proteome_enrichment)
    p_values <- lapply(enrichment_list2, get_p_value)
    p_value_matrix <- make_matrix(disorder_groups = disorder_groups, abundance_groups = abundance_groups)
    p_value_matrix <- fill_matrix(p_values, p_value_matrix)
    p_value_matrix <- p_value_matrix[nrow(p_value_matrix):1 , ]

    effect_size_list <- lapply(enrichment_list, get_effect_size)
    effect_size_matrix <- make_matrix(disorder_groups = disorder_groups, abundance_groups = abundance_groups)
    effect_size_matrix <- fill_matrix(effect_size_list, effect_size_matrix)
    effect_size_matrix  <- effect_size_matrix [nrow(effect_size_matrix ):1 , ]

  return(list(p_value = p_value_matrix,
              effect_size = effect_size_matrix,
              enrichment_list = enrichment_list2,
              binned_genes = binned_genes))
}


stopifnot(FALSE)
master <- read.delim('data/merged_tables/human/actual_master.tsv',
                     sep = "\t")
# log2+1 abundance
master$abundance <- log2(master$abundance+1)

abundance_groups <-  5
disorder_groups <-  5
gradients <-  5
alpha <-  0.01
enrichment_column <- 'is_HI'
disorder_column <- 'aa'
abundance_column <- 'ff'

matrices <- get_matrices(master, abundance_groups = abundance_groups,  disorder_groups = disorder_groups, enrichment_column = enrichment_column)


p_value_matrix <- matrices[[1]]

p_breaks <- make_breaks(unlist(p_value_matrix), alpha = alpha, gradients = gradients)
p_colors <- make_colors(gradients = gradients)
pheatmap(mat = p_value_matrix,
         main = paste0("Abundance:",abundance_groups, " x ", "Disorder:",disorder_groups, " Pvalue heatmap"),
         breaks = c(0, alpha, 1),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = c('red4', 'grey')
)

effect_size_matrix <- matrices[[2]]

pheatmap(mat = effect_size_matrix,
         main = paste0("Abundance:",abundance_groups, " x ", "Disorder:",disorder_groups, " % Enriched Heatmap"),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = brewer.pal(n = 9, name = 'Reds'),
         legend_labels = c("%")
)
