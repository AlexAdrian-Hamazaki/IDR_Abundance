To get the unique_uniprot.txt ids. Ran get_uniprotids.sh. This script will combine the disopred and iupred Ids but as a single list with no duplicates

To get actual uniprot_id mappings, I ran unique_uniprot.txt through https://www.uniprot.org/uniprot/?query=yourlist:M2021102592C7BAECDB1C5C413EE0E0348724B682251F90J&sort=yourlist:M2021102592C7BAECDB1C5C413EE0E0348724B682251F90J&columns=yourlist(M2021102592C7BAECDB1C5C413EE0E0348724B682251F90J),id,entry%20name,reviewed,protein%20names,genes,organism,length
Using the input as "UniProtKB AC/ID" and the output as UniProtKD/Swiss-prot

20252 out of the 20303 mapped to UniProtKD/Swuiss-prot ids! Investingating a couple of the non matched ones looks like they are obsolete Ids

Downloaded the tsv as uniprot_translation.tsv. But I won't actually do anything with that because I plan on joining everything onto these data frames VIA their uniprot ids.
So I don't need the uniprot_symbol info added to this data frame
