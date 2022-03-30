yeast_OE_genes.txt
downloaded from :https://www.sciencedirect.com/science/article/pii/S1097276505018538#app2\
supp table 6 has all toxic genes. Was screened from 85% of the genome
In this table, the third  column represents the growth assessment by spot dilution:

1 = dead
5 = wild-type growth

We don't see any 5 because this table is only for the toxic genes. So we get 4 to 1 scores

yeast_HI_genes.txt
downloaded from: https://academic.oup.com/genetics/article/169/4/1915/6059537?login=true
Supp table 2 and Supp table 1.
I Included both because I want the mutants from supp table 1 as they are still toxic HI phenotypes

awk '{$1=$1}1' yeast_HI* | awk '{print $1"\t"$2"\t"$3}' > yeast_HI_genes.tsv
was used to format the yeast_HI_genes



The Yeast OE and HI Genes are given with their "Systematic Names"

get_unique_ids.sh was run. Giving us unique_systematic_names.txt which is ALL of the systematic names in the Hi and OE dataset

the unique_systematic_names were transformed to Uniprot IDs using https://www.uniprot.org/uploadlists/
Using Ensembl Genome Protein as input and UniprotKB/Swiss-Prot as output

928 Systematic names of the 936 matched to Uniprot Ids

Data was downloaded as "uniprot_translated.tsv"

Data was processed as "entrez_to_uniprot.tsv" using:
awk -F"\t" '{print $1"\t"$2"\t" $6}' uniprot_translated* > entrez_to_uniprot.tsv


new_OE.csv was gotten from sup table 2. Sheet top_HI data. Got uniprot IDs from the
transformer from Ensembl Genomes Proteins Idenifieers to UniprotKB IDs

