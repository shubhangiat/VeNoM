#!/bin/bash

#Sorts matching subgraphs in out_* files according to number of vertices in them and their chi-square values respectively 
#USAGE: ./sort_output.sh filename

lno=`grep -n "The top-[0-9]* subgraphs" $1 |  cut -f1 -d:`
let "lno = $lno + 2"

tail -n +${lno} $1 | sed '/^$/d' | head -n -1 #| awk -F',' 'BEGIN{OFS="\t"} {print NF-1, $0}' | sed 's:==>:==>\t:g' | sort -nt$'\t' -rk1,1 -k3,3 | cut -f2,3 #| cat -n
