#!/bin/bash

awk '{print $1}' *yeast_HI* > systematic_names.txt
awk '{print $1}' *yeast_OE* >> systematic_names.txt

cut -f1 systematic_names.txt | sort | uniq > unique_systematic_names.txt

