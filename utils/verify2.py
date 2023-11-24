# -*- coding: utf-8 -*-
"""

Created on Sun Jun 23 19:56:56 2019

@author: genie

USAGE: python verify2.py query_edge_file target_graph_edge_file sorted_answer_file <n number of results to be verified>

Returns average accuracy of the top-n query graph matches
q = query graph, g = output subgraph, |g| = # edges in g, 
Accuracy for each output subgraph = # edges in g isomorphic to those of q *100/ |q|

"""
import sys
import os
from collections import defaultdict
from numpy import mean

if(len(sys.argv)!=4):
    print(sys.argv)
    print("""\nUSAGE: python verify2.py target_graph_edge_file <file containing paths to query edge file and sorted output file> <n number of results to be verified>

    Returns average accuracy of the top-n query graph matches
    q = query graph, g = output subgraph, |g| = # edges in g, 
    Accuracy for each output subgraph = # edges in g isomorphic to those of q *100/ |q|\n""")
    exit(0)
    
adj = defaultdict(set)  # adjacency list
with open(sys.argv[1]) as f:
    for l in f:
        v1, v2 = l.strip().split()[:2]
        #print(v1 + "\t" + v2)
        if(v1<v2):
            adj[v1].add(v2)
        else:
            adj[v2].add(v1)

#print(len(adj))

with open(sys.argv[2]) as f:
    # read query file and corresponding sorted output file paths
    for fline in f:
        # print(fline)

        # check if file exists
        if not os.path.exists(fline.strip().split()[0]):
            # print(fline.strip().split()[0], "did not exist")
            continue

        # read query edges
        with open(fline.strip().split()[0]) as qry_efile:
            raw = qry_efile.readlines()
        qry_edges = [x.strip().split()[:2] for x in raw]
        # print(qry_edges)

        # check for existence of output file
        if not os.path.exists(fline.strip().split()[1]):
            # print(fline.strip().split()[1], "did not exist")
            continue

        # read top n output subgraph from file
        with open(fline.strip().split()[1]) as sorted_op_file:
            raw = sorted_op_file.readlines()
        raw = raw[:min(len(raw),int(sys.argv[3]))]
        # print(raw)

        # Calculating accuracy
        accuracy = list()
        chisq_list = list()
        while(raw):         # for each output subgraph
            line = raw.pop(0)
            # print(line)
            pairs = line.strip().split(",")[:-1]
            # print(pairs)
            chisq = line.strip().split(",")[-1].replace("==>","").strip()
            # print(chisq)
            # create mapping of target (value) vertex to query (key) vertex
            vmap = dict()
            for p in pairs:
                k = p.split()[1].strip(" ()")
                v = p.split()[0].strip()
                vmap[k] = v
            # check for existence of all query graph edges in output subgraph
            count_iso_edges = 0
            # for each query edge
            for qe in qry_edges:
                # compute mapped vertices
                if(qe[0] not in vmap or qe[1] not in vmap):
                    continue
                v1 = vmap[qe[0]]
                v2 = vmap[qe[1]]
                if(v1 > v2):
                    v1, v2 = v2, v1
                if(v1 in adj and v2 in adj[v1]):
                    count_iso_edges = count_iso_edges + 1
            acc = float(count_iso_edges) / len(qry_edges)
            # print(count_iso_edges)
            # print(len(qry_edges))
            print("Accuracy:\t" + str(acc))
            # print("===========================================================================\n	")
            accuracy.append(acc)
            chisq_list.append(float(chisq))

        max_acc=0.0
        chi_of_max_acc=0.0
        for i in range(len(accuracy)):
            if max_acc<accuracy[i]:
                max_acc = accuracy[i]
                chi_of_max_acc = chisq_list[i]

            
        #print("Average accuracy:\t" + str(round(mean(accuracy), 2)))
        #print("Max accuracy:\t" + str(round(max(accuracy), 2)))
        print(fline.strip().split()[1] + ":\tMax Accuracy:\t" + str(round(max_acc,2)) + "\tChi-sq:\t" + str(round(chi_of_max_acc,4)) + "\tIndex:\t" + str(i))
