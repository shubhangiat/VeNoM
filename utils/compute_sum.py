import sys
from numpy import mean

# Read a comma separated list of values from stdin (pipe the output) and compute it's sum

line=sys.stdin.read()
#print line
if(len(line)==0):
	print(0)
	exit(0)
x=[float(e) for e in line.strip(',').split(",")]
#print(x)
print(len(x), round(sum(x),2))
#print(round(sum(x),2))
