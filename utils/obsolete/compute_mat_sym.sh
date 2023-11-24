#!/usr/bin/env python3
'''
USAGE: ./compute_mat_sym.sh "space separated list of query triplets" "space separated list of input graph triplets"
# triplets MUST not have any space, list must be passed in quotes
'''

import sys


if(len(sys.argv)!=3):
#	print(sys.argv)
	print("USAGE: ./compute_mat_sym.sh \"query triplet\" \"space separated list of input graph triplets\"")
	exit(0)

empty_lab = "phi"
	
q_trip_list = sys.argv[1].strip(", ").split(", ")
ip_trip_list = sys.argv[2].strip(", ").split(", ")

#print(q_trip_list)
#print(ip_trip_list)

#print("\n")
for qt in q_trip_list:
	qt = qt.strip("()")
	q_set = [qt.split(",")[0], qt.split(",")[2]]
	mat_sym_str = ""
	for ipt in ip_trip_list:
#		print(ipt)
		ipt = ipt.strip("()")
		ip_set = [ipt.split(",")[0], ipt.split(",")[2]]
		intersection = []
		for q in q_set:
			if(q in ip_set):
				intersection.append(q)
				ip_set.remove(q)	
		if(empty_lab in intersection):
			intersection.remove(empty_lab)
#		print(str(q_set)+" "+str(ip_set)+" "+str(intersection)+" (s" + str(len(intersection)) + "," +  ipt.split(",")[3] + ") ")
		mat_sym_str = mat_sym_str + "(s" + str(len(intersection)) + "," +  ipt.split(",")[3] + "),"
	print("\"" + mat_sym_str + "\" ")
#print("\n")
