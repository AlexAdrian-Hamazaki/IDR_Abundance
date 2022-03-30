#~~~~~~~~~ Functions ~~~~~~~~~~~~

count_string <- function(vector, string_to_count, type = "") {
  #Count the number of occurances of a string in a data_frame's column
  #Vector is the column name that you want to count in

  # if type == not. Then we count the number of rows that DONT have that string
  if (type != 'not'){
    n <- sum(grepl(pattern = string_to_count, vector))
  }
  if (type == 'not'){
    n <- sum(!grepl(pattern = string_to_count, vector))
  }

  return (n)
}

count_HI <- function(data_frame) {
  #In a dataframe, in the HI column, count how many values are HI, and not HI

  tallied <- count(data_frame, HI)
  HI <- tallied[ (tallied[, 1] == 1 | tallied[, HI] == 2 | tallied[, HI] == 3) , n]
  return(HI)
}

get_highest_percent <- function(data_frame, column_name, high_percent, from = 'high') {
  # Filters a data frame for the highest percent values of a given column name
  # column_name is give as a str
  # percent must be given in proportion format (ex: 0.3 not 30%)

    # if from = low, then you are actually looking at the data from the LOW percent.
  # IF from = high, then you are looking at data from the high percent
  column_name <- enquo(column_name)
  if (from == 'low') {
      data_frame <- arrange(data_frame,(!!column_name))
  } else {
    data_frame <- arrange(data_frame, desc(!!column_name))
  }
  highest_row <- nrow(data_frame)*high_percent

  return(data_frame[1:highest_row, ])

}

get_lowest_percent <- function(data_frame, column_name, high_percent, from = 'high') {
  # Filters a data frame for the highest percent values of a given column name
  # Takes the rows that are NOT in the highest percent
  # Splits those rows in two. The higher will be the middle bin and the lower will be the lowest bin
  # Return: The lowest bin
  # column_name is give as a str
  # percent must be given in proportion format (ex: 0.3 not 30%)


  column_name <- enquo(column_name)
  if (from == 'low') {
      data_frame <- arrange(data_frame,(!!column_name))
  } else {
    data_frame <- arrange(data_frame, desc(!!column_name))
  }
  highest_row <- nrow(data_frame)*high_percent
  mid_low_df <- data_frame[highest_row:nrow(data_frame),]
  middle_index <- nrow(mid_low_df)/2
  low_dataframe <-   mid_low_df[middle_index:nrow(mid_low_df),]

  return(low_dataframe)
}

get_middle_percent <- function(data_frame, column_name, high_percent, from = 'high') {
  # Filters a data frame for the highest percent values of a given column name
  # Takes the rows that are NOT in the highest percent
  # Splits those rows in two. The higher will be the middle bin and the lower will be the lowest bin
  # Return: The lowest bin
  # column_name is give as a str
  # percent must be given in proportion format (ex: 0.3 not 30%)
  column_name <- enquo(column_name)
  if (from == 'low') {
      data_frame <- arrange(data_frame,(!!column_name))
  } else {
    data_frame <- arrange(data_frame, desc(!!column_name))
  }
  highest_row <- nrow(data_frame)*high_percent
  mid_low_df <- data_frame[highest_row:nrow(data_frame),]
  middle_index <- nrow(mid_low_df)/2 + 1
  mid_dataframe <- mid_low_df[1:middle_index,]
  return(mid_dataframe)
}


get_int_counts <- function(data_frame, column_name) {
  # Returns a data frame of 1x3 where you get the total number of rows in data frame,
  # the number of occurances of "1" in a selected column
  # the number of NA values in selected column.

  #column_name must be inputted as a string


  output <- summarize(data_frame,
                      total = n(),
                      counts = sum(!!sym(column_name) , na.rm = TRUE),
                      not_counts = sum(is.na(!!sym(column_name)))
            )
  return(output)
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
  #Purpose: Create an empty list where you will put your genes binned for abundance x disorder
  #Return: Empty list of abundance_groups x disorder_groups length where you should bin your genes

  #@Param:
  #abundance_groups = number of abundance groups (int)
  #disorder_groups = number of disorer groups (int)

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

make_matrix <- function(abundance_groups, disorder_groups) {
  #Purpose: Make a matrix of abundance_groups x disorder_groups dimension
  #Return: Empty matrix of abundance_groups x disorer_groups size

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
create_bins <- function(master, abundance_groups, disorder_groups, abundance_col, disorder_col) {
  # Given a master data frame of genes, bin the genes based on disorder and abundance level
  # Master = Master data frame that contains genes, abundance data, and disorder data
  # abundance_groups = the amount of groups you want your abundance data to be split into
  # disorder_groups = the amount of groups you want your disorder data to be split into
  # disorder_col = the column name of your disorder data
  # abundance_col = the column name of your disorder data
      master <- master %>%
      filter(!is.na(!!(sym(disorder_col)) & !is.na(!!(sym(abundance_col)))))

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