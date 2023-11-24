#!/bin/bash

LC_ALL=C fgrep "Degree computed for $2" output_all
LC_ALL=C fgrep "Chi square was computed for ($2" $1
LC_ALL=C fgrep -cw $2 ~/code/ppi_dataset/protein_edges_perturbed.txt
tail -n +80215 output_all | LC_ALL=C fgrep "ID: $2" | cut -f1,2,4,6,8
