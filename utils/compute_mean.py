import sys
from decimal import *
from numpy import mean

# Read a comma separated list of values from stdin (pipe the output) and compute it's mean

line=sys.stdin.read()
prec = 2
if len(sys.argv)==2:	
	prec = int(sys.argv[1])
#print line
if(len(line)==0):
	print(0)
	exit(0)
x=[Decimal(e) for e in line.strip(',').split(",")]
#print(x)
#print(len(x), round(mean(x),2))
print(round(Decimal(mean(x)),prec))
