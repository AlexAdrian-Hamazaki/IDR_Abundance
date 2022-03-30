#Purpose:
# These functions are required for summarizing and creating a p_value for each of your disorer/abundance x enrichment bins
# After this you are ready for visualization


summarize_enrichment <- function(master, enrichment_column_name) {
  # Purpose: Summarize how many counts of '1' you have in your enerichment_column_name
  # Return: A 1x3 data frame that shows how many '1's(hits) and how many '0's(not_hits) you have in your enricmnent

  #@Param
  #enrichment_column_name = the name of your enrichment column of interest (str): For example, HI/OE or stress_granule enrichment


  enrichment_sum <- summarize(master,
                      total = n(),
                      counts = sum(!!sym(enrichment_column_name), na.rm = TRUE),
                      not_counts = total-counts)

  return(enrichment_sum)
}


get_effect_size <- function(enrichment_summary_table, hit_column = 'counts', miss_column = 'total') {
  # Purpose: Append the effect size onto a enrichment_summary_table
  # For use in an lapply loop
  # Returns: an enrichment summary table with effect size appended

  new_table <- enrichment_summary_table %>%
    mutate(effect_size = !!sym(hit_column)/ !!sym(miss_column))

  return(new_table)
}

get_p_value <- function(enrichment_df, alternative = 'greater', alpha = alpha, hit_column = 'counts', miss_column = 'not_counts') {
  # Purpose: Get P_values from fisher exact test on each of the enrichment summary data frames
  # Returns: List of p_values corresponding to each bin

  # Peforms right tailed fishers exact test.
  # For use in an lapply

  fisher_table <- enrichment_df %>%
    select(!!sym(hit_column), !!sym(miss_column))

  fisher_test <- fisher.test(fisher_table, alternative = alternative, conf.level = 1-alpha)
  p_val <- fisher_test$p.value
  return (p_val)
}

pull_effect_sizes <- function(enrichment_df) {
  #Purpose: Get the effect sizes out of the enrichment_list
  # For use in lapply
  #Returns: The effect size of an enrichment_df
  return(unlist(enrichment_df[ 1 , 4 ]))
}




############################################
#These functions are used exclusively for XX_enrichment_in_bin.R scripts
get_proteome_in_bin <-  function(bin, proteome) {
  # Function: For each bin, identify how many proteins are in that bin. Then use the proteome to identify how many
  # proteins are not in that bin

  # Could be used in lapply

  # Returns: 1x3 data frame summarizing the proteome in a given bin
  bin_size <- nrow(bin)
  proteome_size <- nrow(proteome)
  not_in_bin <- proteome_size - bin_size

  summary_data_frame <- data.frame( total_num = proteome_size,
                                    num_in_bin = bin_size,
                                    num_not_in_bin = not_in_bin)
  return(summary_data_frame)
}

get_enrichment_in_bin <-  function(bin, proteome, enrichment_column_name) {
  # Function: For each bin, identify how many proteins that are of a certain enrichment type.
  # Then use the proteome to identify how many proteins of that enrichment type are not in our bin

  #@PARAM
  # bin = dataframe of genes binned for abundance x disorder
  # proteome = master proteome data frame
  # enrichment_column_name = name of column of whatever enrichment you are looking for (HI or OE or SG, or ...)

  # Example: If enrichment_column_name = HI (happloinsufficiency)
  # In bin, identify how many are HI
  # in proteome identify how many are HI
  # then identify how many are HI and not in our bin

  # Could be used in lapply

  # Returns: 1x3 data frame summarizing the proteome in a given bin


  num_enriched_in_bin <- sum(bin[ ,enrichment_column_name], na.rm = TRUE)
  num_enriched_in_proteome <- sum(proteome[ ,enrichment_column_name], na.rm = TRUE)
  num_not_in_bin <- num_enriched_in_proteome - num_enriched_in_bin

  summary_data_frame <- data.frame( total_num = num_enriched_in_proteome,
                                    num_in_bin = num_enriched_in_bin,
                                    num_not_in_bin = num_not_in_bin)
  return(summary_data_frame)
}

rbind_dfs_in_list <- function(list1, list2) {
  #Function: rbind the first df in list 1 to the first df in list2, then the second, third, etc
  #List1 and list2 must have same length
  #list1 and list2 are lists of data frames

  merged_list <- list()

  for (i in 1:length(list1)) {
    merged_df <- rbind(list1[[i]], list2[[i]])
    merged_list <- append(merged_list, list(merged_df))
  }
  names(merged_list) <- names(list1)
  return(merged_list)
}