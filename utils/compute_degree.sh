#!/usr/bin/env python
'''
USAGE: compute_degree.sh "space separated list of (label,probability)"
# label-probability pair MUST not have anny space, list must be passed in quotes
'''

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 21:00:32 2018

@author: genie
"""

import sys
import itertools

if(len(sys.argv)!=2):
	print("USAGE: compute_degree.sh \"space separated list of (label,probability)\"")
	exit(0)

neigh = sys.argv[1].split()
neigh_pr = []

for n in neigh:
	n = float(n.strip("()").split(',')[1])
	neigh_pr.append(n)

# print(neigh_pr)

#compute power set of indices
pset = set()
deg = [0] * (len(neigh_pr)+1)
for n in range(len(neigh_pr) + 1):
	for sset in itertools.combinations(range(len(neigh_pr)), n):
		cur_deg = len(sset)
#		print cur_deg
		deg_pr = 1
		for i in range(len(neigh_pr)):
			if(i in sset):
				deg_pr = deg_pr*neigh_pr[i]
			else:
				deg_pr = deg_pr*(1-neigh_pr[i])
		deg[cur_deg] = deg[cur_deg] + deg_pr

#		pset.add(sset)
#		print(len(sset))

#deg = [0] * (len(neigh_pr)+1)
#for sset in pset:
#	cur_deg = len(sset)
#	deg_pr = 1
#	for i in range(len(neigh_pr)	):
#		if(i in sset):
#			deg_pr = deg_pr*neigh_pr[i]
#		else:
#			deg_pr = deg_pr*(1-neigh_pr[i])
#	deg[cur_deg] = deg[cur_deg] + deg_pr
#print(deg)

total_deg = 0
for i in range(len(deg)):
	total_deg = total_deg + i*deg[i]
print(round(total_deg,2))
