install.packages('tidyverse')
install.packages('ggplot2')
install.packages('ggdark')
install.packages('RColorBrewer')
install.packages('infer')
install.packages('reshape2')
install.packages('hash')
install.packages('pheatmap')
install.packages('dunn.test')
install.packages('rmarkdown')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")


