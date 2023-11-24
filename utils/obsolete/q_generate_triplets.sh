#!/usr/bin/env python3
'''
USAGE: ./q_generate_triplets.sh v_label "comma separated list of labels"
# list must be passed in quotes
'''

import sys

if(len(sys.argv)!=3):
	print(sys.argv)
	print("\nUSAGE: ./q_generate_triplets.sh v_label \"comma separated list of labels\"")
	exit(0)
	
v_label = sys.argv[1]
neigh = sys.argv[2].split()

triplet_str = ""

while(len(neigh)<2):
	neigh.append("phi")

for i in range(len(neigh)):
	src_lab = neigh[i].strip(',')
	for j in range(i+1, len(neigh)):
		dst_lab = neigh[j].strip(',')
		trip = "(" + src_lab + "," + v_label + "," + dst_lab + ")"
		triplet_str = triplet_str + trip + ", "
	
print(triplet_str)
