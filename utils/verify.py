# -*- coding: utf-8 -*-
"""

Created on Sun Jun 23 19:56:56 2019

@author: genie

USAGE: ./verify.sh query_edge_file target_graph_edge_file sorted_answer_file <n number of results to be verified>

Returns average accuracy of the top-n query graph matches
q = query graph, g = output subgraph, |g| = # edges in g, 
Accuracy for each output subgraph = # edges in g isomorphic to those of q *100/ |q|

"""
import sys
import subprocess as sp
from numpy import mean

if(len(sys.argv)!=5):
    print(sys.argv)
    print("""\nUSAGE: ./verify.sh query_edge_file target_graph_edge_file sorted_answer_file <n number of results to be verified>

    Returns average accuracy of the top-n query graph matches
    q = query graph, g = output subgraph, |g| = # edges in g, 
    Accuracy for each output subgraph = # edges in g isomorphic to those of q *100/ |q|\n""")
    exit(0)
    
# read query edges
with open(sys.argv[1]) as qef:
    raw = qef.readlines()
qe = list()
while(raw):
    qe.append((raw.pop(0)).strip().split())
# print(qe)

# read top n output subgraph from file
with open(sys.argv[3]) as rf:
    raw = rf.readlines()
raw = raw[:min(len(raw),int(sys.argv[4]))]
#print(len(raw)); exit(0)

# Calculating accuracy
accuracy = list()
precision = list()
chisq_list = list()
while(raw):         # for each output subgraph
    line = raw.pop(0)
    # print(line)
    pairs = line.strip().split(",")[:-1]
    chisq = line.strip().split(",")[-1].replace("==>","").strip()
    # create mapping of target vertex to query vertex
    vmap = dict()
    k = map(lambda x: x.split()[1].strip("()"), pairs)  # keys, query vertices
    v = map(lambda x: x.split()[0], pairs)              # values, corresponding target vertex
    for i in range(len(pairs)):
        vmap[k[i]] = v[i]
    # retrieve all edges of output subgraph from target graph
    # for all vertices in v
    nodelist = ""
    for node in v:
        nodelist = nodelist + "-e " + node + " "
    # grep edges from target graph edge file
    grep_op = sp.check_output("ls " + sys.argv[2] + " | xargs -P50 fgrep -w " + nodelist, shell=True)
    tge = map(lambda x: x.split()[:-1], grep_op.splitlines())
    ans_edges = list()
    for e in tge:
        # for vertices in e (for each edge in tge) check complete containment in output subgraph vertices
        # i.e., only consider edges that will belong to the output matching subgraph
        result = all(elem in v for elem in e)
        if result:
            ans_edges.append(e)
#    print(ans_edges)
    # check for existence of all query graph edges in output subgraph
    count_iso_edges = 0
    for kq in qe:    # for each edge in query graph
        # check if all the query vertices have a mapping, skip if not
        if(not all(elem in k for elem in kq)):
            continue
        # find its mapping to target graph
        tv1 = vmap[kq[0]]
        tv2 = vmap[kq[1]]
        # search edge tv1-tv2 in ans_edges
        tve = sorted([tv1, tv2])
        for e in ans_edges:
            if(tve ==  sorted(e)):
                count_iso_edges = count_iso_edges + 1
                break
    acc = float(count_iso_edges) / len(qe)
    prec = 0 if len(ans_edges)==0 else float(count_iso_edges) / len(ans_edges)
#    print(count_iso_edges)
#    print(len(qe))
#    print("Accuracy:\t" + str(acc))
#    print("===========================================================================\n	")
    accuracy.append(acc)
    precision.append(prec)
    chisq_list.append(float(chisq))

max_acc=0.0
prec_of_max_acc=0.0
chi_of_max_acc=0.0
for i in range(len(accuracy)):
    if max_acc<accuracy[i]:
        max_acc = accuracy[i]
        prec_of_max_acc = precision[i]
        chi_of_max_acc = chisq_list[i]

    
#print("Average accuracy:\t" + str(round(mean(accuracy), 2)))
#print("Max accuracy:\t" + str(round(max(accuracy), 2)))
print("Max Accuracy:\t" + str(round(max_acc,2)) + "\tPrecision:\t" + str(round(prec_of_max_acc, 2)) + "\tChi-sq:\t" + str(round(chi_of_max_acc,4)))
for i in range(len(accuracy)):
    print("\t".join(["Accuracy:", str(round(accuracy[i],2)), "Precision:", str(round(precision[i],2)), "Chi-sq:", str(round(chisq_list[i],2))]))

