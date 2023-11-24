#!/usr/bin/env python3
'''
USAGE: ./sum_simlevel.sh <list of "comma-space separated list of (simlevel,probability)">
# simlevel-probability pairs MUST not have any space, lists must be passed in quotes
'''

import sys

if(len(sys.argv)<2):
#	print(sys.argv)
	print("USAGE: ./sum_simlevel.sh <list of \"comma-space separated list of (simlevel,probability)\">")
	exit(0)

#print(sys.argv)
norm_str = ""
observed = [0]*3
prec = 4

for i in range(1,len(sys.argv)):
	pair_list = sys.argv[i].strip(",\"()").split("),(")
#	print(pair_list)
	s = [0] * 3
	for p in pair_list:
#		print(p)
		p_first = p.strip("(), ").split(",")[0]
		p_second = float(p.strip("(), ").split(",")[1])
		if(p_first=="s0"):
			s[0] = s[0]+ p_second
		elif(p_first=="s1"):
			s[1] = s[1]+ p_second
		elif(p_first=="s2"):
			s[2] = s[2]+ p_second
#	print("(" + str(round(s[0],2)) + ", " + str(round(s[1],2)) + ", " + str(round(s[2],2)) + ")")
	sum_s = s[0]+s[1]+s[2]
	norm_str = norm_str + "(" + str(round(s[0]/sum_s,prec)) + ", " + str(round(s[1]/sum_s,prec)) + ", " + str(round(s[2]/sum_s,prec)) + ")\n"
	observed[0] = round(observed[0] +  s[0]/sum_s, prec)
	observed[1] = round(observed[1] +  s[1]/sum_s,prec)
	observed[2] = round(observed[2] +  s[2]/sum_s,prec)
#print("\nNormalized: \n"+ norm_str)
#print("\nObserved: " + str(observed))
print(str(observed).strip("[]"))
#print("\n")
