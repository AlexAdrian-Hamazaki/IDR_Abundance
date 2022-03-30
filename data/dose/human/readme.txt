Census file contails human oncogenes gotten from: 
https://cancer.sanger.ac.uk/census


Clini file contains human HI genes gotten from https://clinicalgenome.org/start/


To get the entrez numbers isolated, and to get the unique ones only. I ran get_unique_ids.sh

These ids put into 
https://www.uniprot.org/uploadlists/

Using "Entrez" as the input language and "UniProtKB/Swiss-prot" as output language

1833 out of 1863 Entrez ids matched to the uniprot Ids. 
This table was downloaded as "entrez_translation.tsv"

Only the two relevant columns were selected. The entrez_id and the UniprotID. This was done via
awk '{print $1"\t"$2}' entrez_translation* > entrez_to_uniprot.tsv


entrez_to_uniprot.tsv  will be used to translate the overexpression and happloinsufficiency data so that they can be merged onto the disorder table
