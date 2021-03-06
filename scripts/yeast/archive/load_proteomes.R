library(Biostrings)
library(utils)
library(stringr)

# Opens up Human and yeast proteome fasta files, wrangles them into a tidy dataframe
# Saves protomes as df_human_proteome.tsv and df_yeast_proteome.tsv

# Open up the Human and Yeast proteomes as FASTA files
yeast_proteome <- Biostrings::readAAStringSet('data/proteomes/uniprot-filtered-reviewed yes+AND+organism Saccharomyces+cerevisiae+(st--.fasta')
human_proteome <- Biostrings::readAAStringSet('data/proteomes/uniprot-proteome UP000005640+reviewed yes.fasta')

# Get Gene names from fasta
sep_name <- function(fasta_rowname) {
  my_split <- str_split(fasta_rowname, pattern = "\\|")
  rough_name <- my_split[[1]][3]
  clean_name <- str_extract(rough_name, pattern = ".*(?=_)")
  return(clean_name)
}
sep_uniprot <- function(fasta_rowname) {
  my_split <- str_split(fasta_rowname, pattern = "\\|")
  rough_name <- my_split[[1]][2]
  return(rough_name)
}
get_seq <- function(fasta_file) {

  l_sequence <- list()

  for (n_row in 1:length(fasta_file)) {
    sequence <-  toString(fasta_file[n_row])
    l_sequence <- append(l_sequence,sequence)
  }
  return(l_sequence)
}

# Transfer Yeast Fasta into a data frame
message("Doing Yeast Proteome")
df_yeast_proteome <- data.frame(matrix(NA, nrow = length(yeast_proteome), ncol = 4))
colnames(df_yeast_proteome) <- c('symbol', 'seq', 'length','uniprot_id')

gene_names <- lapply(names(yeast_proteome), sep_name)
gene_uniprot <- lapply(names(yeast_proteome), sep_uniprot)
sequences <-  get_seq(yeast_proteome)
widths <- width(yeast_proteome)

df_yeast_proteome$symbol <- unlist(gene_names)
df_yeast_proteome$seq <- unlist(sequences)
df_yeast_proteome$length <- unlist(widths)
df_yeast_proteome$uniprot_id <- unlist(gene_uniprot)

# Transfer Human Fasta into a data frame
message("Doing Human Proteome")
df_human_proteome <- data.frame(matrix(NA, nrow = length(human_proteome), ncol = 4))
colnames(df_human_proteome) <- c('symbol','length','uniprot_id', 'seq')

hum_gene_names <- lapply(names(human_proteome), sep_name)
hum_gene_uniprot <- lapply(names(human_proteome), sep_uniprot)
hum_sequences <-  get_seq(human_proteome)
hum_widths <- width(human_proteome)

df_human_proteome$symbol <- unlist(hum_gene_names)
df_human_proteome$seq <- unlist(hum_sequences)
df_human_proteome$length <- unlist(hum_widths)
df_human_proteome$uniprot_id <- unlist(hum_gene_uniprot)

# Write Data Frames
write.table(x = df_yeast_proteome,
            file = "data/proteomes/df_yeast_proteome.tsv",
            quote = FALSE,
            sep = "\t")
write.table(x = df_human_proteome,
            file = "data/proteomes/df_human_proteome.tsv",
            quote = FALSE,
            sep = "\t")