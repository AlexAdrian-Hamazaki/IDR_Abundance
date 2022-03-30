#Purpose:
# These functions are all requred for binning genes by their abundance x disorder values.
# The main executive function of these is create_bins (see the last one in this script)




get_bounds <- function(master,  group_col, group_size){
  #Purpose: Gets the index positions that you will be splitting up your column into

  #group_col is the column name of your target group as a string
  #group_size is the amount of bins you want to bin your column into


  group <- master[ , group_col]
  stopifnot(!any(is.na(group)))
  #If you stop here then there is an na value in your column. These should have been removed preliminarily

  len <- length(group)
  unit <- len/(group_size)

  bounds <- iterate_bounds(group_size, unit)

  return(bounds)
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
add_bin_info <- function(master, bounds, bounds_column, newcolname) {
  #Purpose: Appends a new column to master data frame indicating what bin that gene is in for a given column
  #Returns: A new master data frame that contains info about what bin a protein is for a given metric (abundance/disorder)

  # @PARAM
  # master = master data frame
  # Bounds = the bounds you will split your frame at
  # bounds_column = the column name (string) that you will be using to split via your bounds
  # newcolname = the new column name you will add to indicate which bin a row is in (str)

  # How it works
  # This function takes in your master data frame, arranges it in descending order via a column, splits the df via
  # bounds, and then onto each of the data frames it adds a new column indicating which bin/bounds it is in
  # Then it merges the data frames back together and returns the frame.

  #Splitting Step
  splits <- get_binned_dfs(master, bounds = bounds, bounds_column = bounds_column)

  #Altering Step
  tagged <- add_bin_tag(splits, newcolname = newcolname)

  #Merging Step
  master_new <- bind_rows(tagged)
  return(master_new)
}

add_bin_tag <- function(splits, newcolname) {
  #add a new column and tag every df with "1" in that column
  #newcolname should be the column name you want to add as a str
  for ( i in 1:length(splits) ) {
    splits[[i]] <- mutate(splits[[i]], "{newcolname}" := i)
  }
  return(splits)
}


get_binned_dfs <- function(master, bounds, bounds_column){
  #Orders the data frame by bounds_column, and then splits the data frame into different data frames based on bounds
  #Returns a list of data frames. Each is a different bin.

  #@Param
  #bounds_column = column name that you are splitting as a str
  #bounds = bound indexes you are going to spli at

  master_ordered <- arrange(master, !!sym(bounds_column))
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

get_bins <- function(abundance_groups, disorder_groups) {
  #Purpose: Helper function for create_list(). This actually makes the bins

  bins <- list()
  for ( i in 1:abundance_groups) {
    for ( j in 1:disorder_groups) {
      bins <- append(bins, list(paste0(i,",",j)))
    }
  }
  return(bins)
}

place_rows_in_bins <- function(master, list_of_bins,abundance_bin,disorder_bin) {
  #Purpose: Bin your genes by abundance and disorder
  #Returns: List of data frames where each data frame is a bin of abundance x disorder

  for (row_number in 1:nrow(master)) {
    #message(paste('Evaluating row number:', row_number))
    row <- master[row_number,]
    binval <- paste0(row[,abundance_bin],",",row[disorder_bin])
    for ( listname in names(list_of_bins)) {
        if ( binval == listname & !is.na(list_of_bins[binval]) ) {
          list_of_bins[[binval]]<- rbind(list_of_bins[[binval]], row)
          break
        } else if ( binval == listname & is.na(list_of_bins[binval]) ) {
          list_of_bins[binval]<- list(row)
          break
        }
    }

  }
  return (list_of_bins)
}

create_bins <- function(master, abundance_groups, disorder_groups, abundance_col, disorder_col) {
  # Given a master data frame of genes, bin the genes based on disorder and abundance level
  # Returns: List of data frames corresponding to the different bins. The x value corresponds to abundance bin. The y value coresponds to disorder bin.


  # Master = Master data frame that contains genes, abundance data, and disorder data
  # abundance_groups = the amount of groups you want your abundance data to be split into
  # disorder_groups = the amount of groups you want your disorder data to be split into
  # disorder_col = the column name of your disorder data
  # abundance_col = the column name of your disorder data
  master <- master %>%
    filter(!is.na(!!sym(abundance_col))) %>%
    filter(!is.na(!!sym(disorder_col)))

    abundance_bounds <- get_bounds(master, group_col = abundance_col, group_size = abundance_groups)
    disorder_bounds <- get_bounds(master, group_col = disorder_col, group_size = disorder_groups)

  #Add a new column indicating what abundance bin each gene is in
  master <- add_bin_info(master,
                           bounds = abundance_bounds,
                           bounds_column = abundance_col,
                           newcolname = 'abundance_bin')
  #Add a new column indicating what disorder bin each gene is in
  master <- add_bin_info(master,
                        bounds = disorder_bounds,
                        bounds_column = disorder_col,
                        newcolname = 'disorder_bin')

  #Create empty list where you will put your bined genes
  list_of_bins <- create_list(abundance_groups, disorder_groups)

  #Place genes into their respective disorder x abundance bin
  binned_genes <- place_rows_in_bins(master, list_of_bins, 'abundance_bin', 'disorder_bin')
  return(binned_genes)
}

