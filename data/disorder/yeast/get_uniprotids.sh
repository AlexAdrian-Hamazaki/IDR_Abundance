#!/bin/bash

awk '{print $1}' *yeast_dis* > uniprots.txt
awk '{print $1}' *yeast_iu* >> uniprots.txt

cut -f1 uniprots.txt | sort | uniq > unique_uniprots.txt

