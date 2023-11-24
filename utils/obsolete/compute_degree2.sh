#!/bin/bash

#USAGE: ./compute_degree2.sh space separated list of probabilities of neighbouring vertices

IFS=' ' read -r -a neigh_pr <<< "$@"
len_npr=`wc -w <<< $@`

pow_set_size=`echo "2^$len_npr" | bc -l`

#echo $pow_set_size
degree=0
for i in `seq 0 $((pow_set_size-1))`
do
	cur_deg=0
	cur_deg_pr=1
	chk_str=""
	for j in `seq 0 $((len_npr-1))`
	do
		pow=`echo "2^$j" | bc -l`
#		echo $i $pow $((i & pow))
		if [ $((i & pow)) -ne 0 ]
		then
#			echo ${neigh_pr[j]}
			chk_str=$chk_str" "${neigh_pr[j]}
			cur_deg=$((cur_deg+1))
			cur_deg_pr=`echo "$cur_deg_pr*${neigh_pr[j]}" | bc -l`
		else
			rev_npr=`echo "1-${neigh_pr[j]}" | bc -l`
#			echo $rev_npr
			chk_str=$chk_str" "$rev_npr
			cur_deg_pr=`echo "$cur_deg_pr*$rev_npr" | bc -l`
		fi
	done
#	echo "$cur_deg :$chk_str"
	degree=`echo "$degree+$cur_deg*$cur_deg_pr" | bc -l`
done

echo $degree
