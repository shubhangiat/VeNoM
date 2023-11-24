import sys
from numpy import mean

# Read a comma separated list of values from stdin (pipe the output) and compute it's MRR for a list length of 20
# add 0s for missing values

def mean_MRR(x, len_y):
    x = [e if e<=10 else 0 for e in x]
    y = [1.0/e if e!= 0 else 0 for e in x] 
    #len_y = input("\nTotal number of observations to compute mean from?\t")
    while(len(y)<len_y):
        y.append(0)
    #print((len(y), mean(y)))
    return mean(y)

line=sys.stdin.read()
#print line
if(len(line)==0):
    print(0)
    exit(0)
x=[float(e) for e in line.strip(',').split(",")]
#print(x)
#print(len(x), round(mean(x),2))
print(round(mean_MRR(x, int(sys.argv[1])),2))
