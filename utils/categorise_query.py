#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 16:45:30 2022

@author: Shubhangi Agarwal

Description:
Categorize a given list of queries as dense or sparse
It is better to supply path of exact query files only
If the query has edges such that the avg. degree of the graph is >=3 => dense

USAGE: 
python categorise_query.py <arg1> <arg2> <arg3>
arg1 - path to file containing list of queries
arg2 - output filename for dense query
arg3 - output filename of sparse query
Format of input file - "absolute path of query vertex file" "absolute path of query edge file"
Format of output file - "path of sorted output file", e.g. - "q_v5/sorted_output_v5_1"
Output format of arg2 and arg3 - "path of query label file" "path of query edge file"
"""

###############################################################################

import sys
import networkx as nx

# CONSTANTS
threshold_degree = 3

if len(sys.argv) != 4:
    print("""Description:
        Categorize a given list of queries as dense or sparse
        It is better to supply path of exact query files only
        If the query has edges such that the avg. degree of the graph is >=3 => dense
    """)
    print("""USAGE: 
        python categorise_query.py <arg1> <arg2> <arg3>
        arg1 - path to file containing list of queries
        arg2 - output filename for dense query
        arg3 - output filename of sparse query
        Format of input file - "absolute path of query vertex file" "absolute path of query edge file"
        Format of output file - "path of sorted output file", e.g. - "exact/q_v5/sorted_output_v5_1"
        Output format of arg2 and arg3 - "path of query label file" "path of query edge file"
    """)
    exit(0)

with open(sys.argv[1]) as f:
    raw = f.readlines()
    inp_data = [x.strip() for x in raw]

qdense = []
qsparse = []

for r in inp_data:
    fqlabel = r.split()[0]
    fqedge = r.split()[1]

    with open(fqedge) as f:
        raw = f.readlines()
    
    qedges = [x.strip().split() for x in raw]

    qg = nx.Graph()
    qg.add_edges_from(qedges)
    n = qg.number_of_nodes()
    e = qg.number_of_edges()

    max_edges = n*(n-1)/2
    threshold_edges = n*1.0*threshold_degree/2       # ensuring average degree of graph is >=3

    # path of sorted output file
    sof_path = fqlabel.split("/")[-2] + "/" + fqlabel.split("/")[-1].replace("q", "sorted_output").replace("label_", "")

    if e >= threshold_edges:
        qdense.append(sof_path)
    else:
        qsparse.append(sof_path)

with open(sys.argv[2], "w") as fdense:
    fdense.write("\n".join(qdense))

with open(sys.argv[3], "w") as fsparse:
    fsparse.write("\n".join(qsparse))
