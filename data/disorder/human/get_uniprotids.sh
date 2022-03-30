#!/bin/bash

awk '{print $1}' human_dis* > uniprots.txt
awk '{print $1}' human_iu* >> uniprots.txt

cut -f1 uniprots* | sort | uniq > unique_uniprots.txt

