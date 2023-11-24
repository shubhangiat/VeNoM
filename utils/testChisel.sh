#!/bin/bash

qry_edge=dataset/q_v4_edge2.txt
./subgraph dataset/ip_v27_label.txt dataset/ip_v27_edgePr.txt <(echo dataset/q_v4_label.txt $qry_edge) > tmp

~/sm_frmwrk/utils_chisel/sort_output.sh tmp > test

echo "Input graph: `head -n1 tmp | sed 's:^.*(::; s:).*$::'`"
echo "Query graph: `/bin/grep 'Query files' tmp | sed 's:^.*\t\(.*\)\t:\1, :'`"

python ~/sm_frmwrk/utils_chisel/verify3.py dataset/ip_v27_edgePr.txt <(echo $qry_edge test ) 10

cat test
