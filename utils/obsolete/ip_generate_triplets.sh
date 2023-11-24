#!/usr/bin/env python3
'''
USAGE: ./ip_generate_triplets.sh v_label "space separated list of (label,probability)"
# label-probability pairs MUST not have any space, list must be passed in quotes
'''

import sys	

if(len(sys.argv)!=3):
#	print(sys.argv)
	print("USAGE: ./ip_generate_triplets.sh v_label \"space separated list of (label,probability)\"")
	exit(0)

v_label = sys.argv[1]
neigh = sys.argv[2].split()

triplet_str = ""

prec = 4

# specified labels exist as neighbours
for i in range(len(neigh)):
	src_lab = neigh[i].strip("()").split(',')[0]
	src_pr = float(neigh[i].strip("()").split(',')[1])
	for j in range(i+1, len(neigh)):
		dst_lab = neigh[j].strip("()").split(',')[0]
		dst_pr = float(neigh[j].strip("()").split(',')[1])
		trip_pr = round(src_pr * dst_pr, prec)
		trip = "(" + src_lab + "," + v_label + "," + dst_lab + "," + str(trip_pr) + ")"
		triplet_str = triplet_str + trip + ", "
		
	# only one label exist as neighbour
	trip_pr = src_pr
	for j in range(len(neigh)):
		if(j!=i):
			trip_pr = trip_pr * (1 - float(neigh[j].strip("()").split(',')[1]))
	trip_pr = round(trip_pr, prec)
	if(trip_pr > 0):
		dst_lab = "phi"
		trip = "(" + src_lab + "," + v_label + "," + dst_lab + "," + str(trip_pr) + ")"
		triplet_str = triplet_str + trip + ", "

# when no neighbour label exists
trip_pr = 1
for j in range(len(neigh)):
	trip_pr = trip_pr * (1 - float(neigh[j].strip("()").split(',')[1]))
trip_pr = round(trip_pr, prec)
if(trip_pr > 0):
	src_lab = "phi"
	dst_lab = "phi"
	trip = "(" + src_lab + "," + v_label + "," + dst_lab + "," + str(trip_pr) + ")"
	triplet_str = triplet_str + trip
print(triplet_str)
		
