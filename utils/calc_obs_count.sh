#!/bin/bash

#USAGE: ./calc_obs_count.sh vertex_pair_label "space separated neighbour list of (label,probability) for query vertex" "space separated neighbour list of (label,probability) for input graph vertex"
# label-probability pairs MUST not have any space, list must be passed in quotes


v_label=$1
q_neigh=$2
ip_neigh=$3

#echo $v_label
#echo $q_neigh
#echo $ip_neigh

# generating triplets
q_trip=`./q_generate_triplets.sh $v_label "$q_neigh"`
ip_trip=`./ip_generate_triplets.sh $v_label "$ip_neigh"`

#echo $q_trip
#echo $ip_trip

# computing matching symbols
mat_sym=`./compute_mat_sym.sh "$q_trip" "$ip_trip"| tr '\n' ' '`

#echo $mat_sym
#eval 'for m in '$mat_sym'; do echo \"$m\"; done'


# computing observed count
obs=`./sum_simlevel.sh $(for m in $mat_sym; do echo $m; done)`

echo -n $obs
