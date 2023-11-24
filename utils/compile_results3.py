#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 10:52:11 2019

@author: Shubhangi

USAGE: ./compile_results <dir path where all the result files are stored>
Compiles query processing time and heap size data
On screen output

"""

import sys
import os
from numpy import mean 

if(len(sys.argv)!=2):
        print("USAGE: ./compile_results2.py <dir path where all the result files are stored>")
        print("Compiles query processing time and heap size data")
        print("On screen output")
        exit(0)
        
dir_path = sys.argv[1]

ph_size = list()        # Primary Heap size
sh_size_max = list()    # Max Secondary Heap size
sh_size_avg = list()    # Avg Secondary Heap size
time_taken = list()     # Time taken
index_time = list()     # If included in output file being read
for f in os.listdir(dir_path):
    #print(dir_path + "/" + f)
    if("sorted" in f or "verify" in f):
        continue
    if(f.endswith(("_1", "_2", "4", "5", "7", "9", "_10", "_14", "_15", "_16", "_20"))):
        continue
    else:
        rf = open(dir_path + "/" + f)
        raw = rf.readlines()
        ph_s = int (raw[7].strip().split(" = ")[1])
        if("NO MATCHING SUBGRAPH FOUND!!!" not in raw[7]):
            # print(raw[8])
            sh_sm = int (raw[8].strip().split(" = ")[1])
            sh_sa = float (raw[9].strip().split(" = ")[1])
            #print(raw[1] + raw[2])
            tt = float (raw[11].strip().split()[-2]) + float (raw[2].strip().split()[-2])   # Query processing time + Time taken to compute subgraphs
        else:
            tt = float (raw[11].strip().split()[-2]) + float (raw[2].strip().split()[-2])   # Query processing time + Time taken to compute subgraphs
        ph_size.append(ph_s)
        sh_size_max.append(sh_sm)
        sh_size_avg.append(sh_sa)
        time_taken.append(tt)
        rf.close()
ph_size.sort()
sh_size_max.sort()
sh_size_avg.sort()
time_taken.sort()
print(dir_path.split("/")[-1])
print("Avg Primary-Heap-size:\t" + str(round(mean(ph_size), 2)))
print("Avg Max-Secondary-Heap-size:\t" + str(round(mean(sh_size_max), 2)))
print("Avg Avg-Secondary-Heap-size:\t" + str(round(mean(sh_size_avg), 2)))
print("Avg time taken:\t" + str(round(mean(time_taken), 2)) + "\tsec")
if(len(time_taken)>9):
    trim = int(0.1*len(time_taken))
    print("\nTrimmed")
    print("Avg Primary-Heap-size:\t" + str(round(mean(ph_size[trim:-1*trim]), 2)))
    print("Avg Max-Secondary-Heap-size:\t" + str(round(mean(sh_size_max[trim:-1*trim]), 2)))
    print("Avg Avg-Secondary-Heap-size:\t" + str(round(mean(sh_size_avg[trim:-1*trim]), 2)))
    print("Avg time taken:\t" + str(round(mean(time_taken[trim:-1*trim]), 2)) + "\tsec")
