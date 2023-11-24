#!/usr/bin/env python
###########
'''
Divide "output_all" file by the number of query graph searches it has
Call from the folder where you want the answer files
'''
import sys
import subprocess as sp

if(len(sys.argv)!=3):
	print("USAGE: " + sys.argv[0] + " <filename to be divided> <number of files it needs to be divided into>")
	exit(0)

# Calling grep for line numbers and corresponding output filenames

n = int(sys.argv[2]) *2
#print(n)
ln_gr_output = sp.Popen(['LC_ALL=C fgrep -nm' + str(n) + ' -e "Query files:" -e "Input graph:" ' + sys.argv[1]], shell=True, stdout=sp.PIPE, stderr=sp.STDOUT)

output_lines = list()
for l in ln_gr_output.stdout.readlines():
#	print(l.strip())
	output_lines.append(l.strip())

for i in range(len(output_lines)-1):
	l1 = output_lines[i]
	start = int(l1.split(":")[0])
	l2 = output_lines[i+1]
	end = int(l2.split(":")[0]) - 1
	ofile_name = "output_" + l1.split("/")[-1].strip(".txt")
#	print(start, end, ofile_name)
	file_output = sp.Popen(["sed -n '" + str(start) + ", " + str(end) + "p' " + sys.argv[1] + " > " + ofile_name], shell=True, stdout=sp.PIPE, stderr=sp.STDOUT)
