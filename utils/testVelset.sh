#!/bin/bash

./subgraph dataset/ip_v27_label.txt dataset/ip_v27_edge.txt <(echo dataset/q_v4_label.txt dataset/q_v4_edge2.txt ) > tmp

~/sm_frmwrk/utils/sort_output.sh tmp > test

echo "Input graph: `head -n1 tmp | sed 's:^.*(::; s:).*$::'`"
echo "Query graph: `/bin/grep 'Query files' tmp | sed 's:^.*\t\(.*\)\t:\1, :'`"

python ~/sm_frmwrk/utils/verify2.py dataset/ip_v27_edge.txt <(echo dataset/q_v4_edge2.txt test ) 10

cat test