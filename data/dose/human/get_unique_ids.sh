#!/bin/bash

awk -F"\t" '{print $3}' Cens** > entrez.txt
awk -F"\t" '{print $2}' Clin* >> entrez.txt

cut -f1 entrez.txt* | sort | uniq > unique_entrez.txt

