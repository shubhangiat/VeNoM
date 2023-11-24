#!/usr/bin/env python

# USAGE: ./compute_degree_pr_lx.sh "space separated neighbour list of (label,probability) for input graph vertex"
# Compute degree of the vertex and different probabilities of existence (of instances) for neighbour labels
# Probability that no instance of a label exists, exactly one instance exists and at least one instance exists

import sys

if(len(sys.argv)!=2):
    #print(sys.argv)
    print("USAGE: ./compute_degree_pr_lx.sh \"space separated neighbour list of (label,probability) for input graph vertex\"")
    exit(0)

input = sys.argv[1]

degree = 0
uniq_lab_pr_inst = dict()
pr_lx_0 = dict()
pr_lx_1 = dict()
pr_lx_2 = dict()

# check if input is comma+space separated
if("), " in input):
    input=input.replace("), ", ") ")
    #print(input)

for pair in input.split():
    #print(pair)
    label = pair.split(',')[0].strip("(")
    pr = float(pair.split(',')[1].strip(")"))
    #print(label + " " + pr)
    degree = degree + pr
    if(label in uniq_lab_pr_inst.keys()):
        l = uniq_lab_pr_inst[label]
        l.append(pr)
        pr_0 = pr_lx_0[label]
        pr_0 = pr_0 * (1-pr)
    else:
        l = [pr]
        pr_0 = (1-pr)
    uniq_lab_pr_inst[label] = l
    pr_lx_0[label] = pr_0
    #print(label + ": ", end=' ')
    #print(uniq_lab_pr_inst[label])
print(round(degree,4))

for label in uniq_lab_pr_inst.keys():
    pr_1 = 0
    if(pr_lx_0[label]!=0):
        for e in uniq_lab_pr_inst[label]:
            pr_1 = pr_1 + (pr_lx_0[label]*e/(1-e))
    else:
        for i in range(len(uniq_lab_pr_inst[label])):
            pr_2 = 1
            for j in range(len(uniq_lab_pr_inst[label])):
                if(i!=j):
                    pr_2 = pr_2 * (1-uniq_lab_pr_inst[label][j])
            pr_1 = pr_1 + (pr_2 * uniq_lab_pr_inst[label][i])
    pr_lx_1[label] = pr_1
    pr_lx_2[label] = 1 - pr_lx_0[label]

for k in pr_lx_0.keys():
    print("(" + k + "," + str(round(pr_lx_0[k], 4)), end=") ")
print()
for k in pr_lx_1.keys():
    print("(" + k + "," + str(round(pr_lx_1[k], 4)), end=") ")
print()
for k in pr_lx_2.keys():
    print("(" + k + "," + str(round(pr_lx_2[k], 4)), end=") ")
print()