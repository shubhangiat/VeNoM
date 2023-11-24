# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 16:33:55 2019

@author: Shubhangi

USAGE: python split_ouput_all.py <output file>

Call from folder where the separated files are to be stored, with the folder structure already in place
The output file contains output over queries which are run one at a time
(target graph is read as many times as the number of query graph files)
The splitting of outputs is done after the line 'Input Graph:' is read
Assumption: The last line of program is 'Input Graph:', the next line is the start of output for the next query graph
Output file naming is done based on edge file name in the output file
e.g. output for query edge file q_vXX_edges_YY is stored in folder q_vXX under the file name output_vXX_YY
Existing files will be overwritten
"""

import sys
from pathlib import Path

if len(sys.argv)!=2:
    print("""\nUSAGE: python split_ouput_all.py <output file>

    Call from folder where the separated files are to be stored, with the folder structure already in place
    The output file contains output over queries which are run one at a time
    (target graph is read as many times as the number of query graph files)
    The splitting of outputs is done after the line 'Input Graph:' is read
    Assumption: The last line of program is 'Input Graph:', the next line is the start of output for the next query graph
    Output file naming is done based on edge file name in the output file
    e.g. output for query edge file q_vXX_edges_YY is stored in folder q_vXX under the file name output_vXX_YY
    Existing files will be overwritten\n""")
    exit(0)
    
with open(sys.argv[1]) as f:
    raw = f.readlines()

while ("Query files:" not in raw[0]):
    raw.pop(0)

while(raw):
    # print(raw[0].strip())
    # Comprehend output file path
    efname = raw[0].strip().split("/")[-1]
    # print(efname)
    dir_name = raw[0].strip().split("/")[-2]
    # print(dir_name)
    of_name = efname.split("/")[-1].replace("q", "output").replace("_edges","")
    # print(of_name)
    parent_dir= raw[0].strip().split("/")[-3]
    # print(parent_dir)
    o_path = f'{parent_dir}/{dir_name}'
    Path(o_path).mkdir(parents=True, exist_ok=True)
    with open(o_path+"/"+of_name, "w") as f:  
        line = raw.pop(0)
        while("Input graph" not in line and raw):
            f.write(line)
            line = raw.pop(0)
        f.write(line)
    while(raw and raw[0].strip()==""):
        raw.pop(0)
