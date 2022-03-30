# Purpose: Problem is that the uniprot proteome has protein symbols and codes
# The abundance data uses string_ids to link to proteins. But the protein symbols and codes are not the same

# To solve this issue, We will be uses a file downloaded from stringdb that has all known alises for a string db id.
# We will then create a df that contains all the string dbs codes , and known affiliated protein symbols

# Output: data/abundance/df_pro_linker.tsv is able to link both the uniprot and stringdb protein names.

#this will be done with protein.alises downloaded from STRING:
#https://string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Homo+sapiens

#NON-FUNCTIONAL 2021-09-13

library(stringr)
library(dplyr)
library(tidyr)
library(rlang)

message('Opening Large Files')
string_db <- read.delim(file = "data/abundance/9606.protein.info.v11.5.txt",
                           header = TRUE,
                           sep = "\t",
                           quote = "",
                           fill = FALSE)

uniprot_db <- read.delim("data/proteomes/df_human_proteome.tsv")
aliases <- read.delim(file = "data/abundance/9606.protein.aliases.v11.5.txt")

get_correct_alias <- function(df_aliases, uniprot_ids) {
  # for a given string db code in the alises database, regex to see if the uniprot or string protein code aliases
  # are affiliated with the string code.

  #output a table that links the uniprot pc codes to string pc codes to a string db code
  message('')
  message(paste(df_aliases[1,1], 'in progress...'))
  output_list <- list(NA,NA)
  names(output_list) <- c('string_id', 'uniprot_id')
  output_list[['string_id']] <- df_aliases[[1,1]]

  tryCatch(
  {#Trycatch open. Should have all for loop

  for (alias in df_aliases$alias) {
    if (any(str_detect(string = uniprot_ids, pattern = alias))) {
      output_list[['uniprot_id']] <- append(output_list[['uniprot_id']], alias)
    }
  } # For loop end
    if (is.na(output_list[['uniprot_id']][1])) {
  output_list [['uniprot_id']] <- output_list[['uniprot_id']][-1]
    }
  }
    ,
  error = function(cond) {
      message(paste("Terminal Error in:", alias))
      message(cond)
    }
  )#End try catch
  return(output_list)
}

filter_out_obv_negs <- function(df_from_l_aliases){
  #For each df in l_aliases, we will remove rows with alises which we do not expect to be protein ids.

  processed <- df_from_l_aliases %>%
    distinct(alias, .keep_all = TRUE) %>%
    #all uniprot ids are either 6 or 10 of length
    filter(str_length(alias) == 6 | str_length(alias) == 10)
  if(is_empty(processed)) {
    processed <- data.frame(status = "EMPTY")
  }
  return(processed)
}

empty_to_NA<- function(vector_length1) {
  if (length(vector_length1)==0) {
    return(NA)
  } else {
    return(vector_length1)
  }
}

omit_wrap <- function(df_aliases, uniprot_ids) {
  if (nrow(df_aliases) == 0) {
    warning(paste(df_aliases[1,1], 'Failed for an unknown reason'))
    return (data.frame(string_id <- df_aliases[1,1],
                       uniprot_id <- NA))
  } else {
    return (get_correct_alias(df_aliases, uniprot_ids))
  }
}
# Group alises by their string id. And seperate them into dfs
l_aliases <- aliases %>%
  group_by(X.string_protein_id) %>%
  group_split()
l_names <- aliases %>%
  distinct(X.string_protein_id)
names(l_aliases) <- unlist(l_names)



#process l_alises to remove false positives, this saves time in the regex later on
message("Removing aliases from our alias file that are almost certainly false positives")
if (!file.exists('../../data/abundance/human/human_aliases.RData')) {
  pro_l_aliases <- lapply(X =  l_aliases,
                          filter_out_obv_negs)
  save(pro_l_aliases, file = '../../data/abundance/human/human_aliases.RData')
} else {
  load('../../data/abundance/human/human_aliases.RData')
}


# Get a linker table
message('Grepping through alises')
raw_linker <- lapply (X = pro_l_aliases,
                  FUN = omit_wrap ,
                  uniprot_ids = uniprot_db$uniprot_id)
message('All Aliases succesfully collapsed')

indexes <- c()
for (n_list in 1:length(raw_linker)) {
  if (length(raw_linker[[n_list]]) != 2) {
    indexes <- append(indexes, n_list)
}}
pro_linker <- raw_linker[-c(indexes)]




df_raw_linker <- as.data.frame(do.call(rbind, pro_linker))
df_raw_linker$string_id <- unlist(df_raw_linker$string_id)

string_db <- rename(string_db, 'string_id' = 'X.string_protein_id')
string_db <- rename(string_db, 'stringdb_symbol' = 'preferred_name')


#Join by string_id
df_raw_linker <- left_join(string_db, df_raw_linker, by = "string_id")
df_raw_linker$uniprot_id <- lapply(df_raw_linker$uniprot_id, empty_to_NA)

#Then join uniprot info on
uniprot_db <- rename(uniprot_db, 'uniprot_symbol' = 'symbol')
uniprot_db$uniprot_id <- as.list(uniprot_db$uniprot_id)

df_raw_linker <- left_join(uniprot_db, df_raw_linker, by = "uniprot_id")

df_pro_linker <- df_raw_linker[,c(5,6,3,1,7,4)]
df_pro_linker <- apply(df_pro_linker, 2 , as.character)
df_pro_linker <- df_pro_linker %>%
  distinct(uniprot_symbol, .keep_all = TRUE)

#Most of the genes in df_pro_bad look like predicted genes and such. I'm ok with them not being here. It looks like
# for the most part, there are no string_ids. which makes sense because they prob wouldn't include them
# in the string functional analysis and stuff. So I'll be happy with my get_correct_alias function for now


#
# #Make columns characters not lists so we can write
message('Saved pro_linker.tsv in data/linkers')
dir.create('data/linkers', showWarnings = FALSE)
write.table(x = df_pro_linker, file = "data/linkers/human_pro_linker.tsv",
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)
#
# df_raw_linker <- apply(df_raw_linker, MARGIN = 2, as.character)
#
# write.table(x = df_raw_linker, file = "data/abundance/human_raw_linker.tsv",
#             col.names = TRUE,
#             sep = "\t",
#             quote = FALSE)


