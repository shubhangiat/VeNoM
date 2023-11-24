# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 11:56:12 2019

@author: genie

Outputs show_name-person_name mapping to screen (redirect output to file)
Edges corresponding to person_id (nconst) not present in name_basics are ignored
Use the three files to get who is connected to which films/ shows
All the files are assumed to be tab separated and have one line headers

"""

import os
import sys

if(len(sys.argv)!=4):
    print("USAGE: creates_dataset.sh title_akas title_prinicpals name_basics")
    exit(0)
    
title_akas = open(sys.argv[1], 'r')
title_principals = open(sys.argv[2], 'r')
name_basics = open(sys.argv[3], 'r')

#############
# Reading input files
#############

show_title=dict()

header = True
for l in title_akas:
    if(header):
        header = False
        continue
#    print(l)
    title_id = l.strip().split("\t")[0]
    title = l.strip().split("\t")[1]
    show_title[title_id] = title

cast=dict()

header = True
for l in title_principals:
    if(header):
        header = False
        continue
#    print(l)
    title_id = l.strip().split("\t")[0]
    person_id = l.strip().split("\t")[1]
    if(title_id not in show_title.keys()):      # Discarding mapping for titles not in title_akas
        continue
    if(title_id not in  cast.keys()):
        cast[title_id] = []
    cast[title_id].append(person_id)
    

person=dict()

header=True
for l in name_basics:
    if(header):
        header = False
        continue
#    print(l)
    person_id = l.strip().split("\t")[0]
    name = l.strip().split("\t")[1]
    person[person_id] = name
    
    
#############
# 'Joining' the three dictionaries to map title of film/show and name of related people
#############
    
for title_id in cast.keys():
    title = show_title[title_id]
    for person_id in cast[title_id]:
        if(person_id in person.keys()):         # Ignore if person_id to name mapping is not present
            print(title + "," + person[person_id])
    