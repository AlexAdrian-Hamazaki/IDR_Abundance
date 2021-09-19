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

hum_with_symbol <- read.delim(file = "data/abundance/hum_with_symbol.tsv")
df_human_proteome <- read.delim(file = "data/proteomes/df_human_proteome.tsv")
aliases <- read.delim(file = "data/abundance/9606.protein.aliases.v11.5.txt")


# Group alises by their string id. And seperate them into dfs
l_string_ids <- aliases %>%
  group_by(X.string_protein_id) %>%
  group_split()
l_names <- aliases %>%
  distinct(X.string_protein_id)
names(l_string_ids) <- unlist(l_names)

get_correct_alias <- function(df_aliases, v_uniprot_aliases, v_string_aliases) {
  # for a given string db code in the alises database, regex to see if the uniprot or string protein code aliases
  # are affiliated with the string code.

  #output a table that links the uniprot pc codes to string pc codes to a string db code
  message('')
  message(paste(df_aliases[1,1], 'in progress...'))
  output_list <- list(NA, NA, NA)
  names(output_list) <- c("string_id", "string_symbol", "uniprot_symbol")

  output_list[['string_id']] <- df_aliases[[1,1]]

  unique_alises <- unique(df_aliases$alias)

  tryCatch(
  {#Trycatch open. Should have all for loop

  for (alias in unique_alises) {
    alias <- str_remove(alias, pattern = "_human")
    if (any(str_detect(string = v_string_aliases, pattern = alias))) {
      output_list[['string_id_alias']] <- alias
    }
    if (any(str_detect(string = v_uniprot_aliases, pattern = alias))) {
      output_list[['uniprot_alias']] <- alias
    }
  } # For loop end

  },
  error = function(cond) {
      message(paste("Terminal Error in:", alias))
    }
  )#End try catch
  return(output_list)
}


raw_linker <- lapply (X = l_string_ids,
                  FUN = get_correct_alias ,
                  v_uniprot_aliases = df_human_proteome$symbol,
                  v_string_aliases = hum_with_symbol$symbol)
df_raw_linker <- as.data.frame(do.call(rbind, raw_linker))

df_pro_linker <- df_raw_linker %>%
  filter(!is.na(string_id_alias) & !is.na(uniprot_alias))

#Most of the genes in df_pro_bad look like predicted genes and such. I'm ok with them not being here. It looks like
# for the most part, there are no string_ids. which makes sense because they prob wouldn't include them
# in the string functional analysis and stuff. So I'll be happy with my get_correct_alias function for now
df_pro_bad<- df_raw_linker %>%
  filter(is.na(string_id_alias) | is.na(uniprot_alias))

print(paste(nrow(df_pro_linker), 'symbols are linked'))
print(paste(nrow(df_pro_bad), 'symbols are not linked'))


#Make columns characters not lists so we can write
df_pro_linker <- apply(df_pro_linker, 2, as.character)
write.table(x = df_pro_linker, file = "data/abundance/human_pro_linker.tsv",
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)

df_raw_linker <- apply(df_raw_linker, MARGIN = 2, as.character)

write.table(x = df_raw_linker, file = "data/abundance/human_raw_linker.tsv",
            col.names = TRUE,
            sep = "\t",
            quote = FALSE)


