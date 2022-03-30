To get the unique_uniprot.txt ids. Ran get_uniprotids.sh.
This script will combine the disopred and iupred Ids but as a single list with no duplicates

To get actual uniprot_id mappings,
I ran unique_uniprot.txt through https://www.uniprot.org/uniprot/
Using the input as "UniProtKB AC/ID" and the output as UniProtKD/Swiss-prot

6729 out of the 6721 mapped. Meaning there were some duplicates. It looks like the duplicates are mostly uncharacterized proteins
that likely have unstable names

Downloaded the tsv as uniprot_translation.tsv.
But I won't actually do anything with that because I plan on joining everything onto these data frames VIA their uniprot ids.
So I don't need the uniprot_symbol info added to this data frame
