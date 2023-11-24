#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 00:33:40 2022

@author: Shubhangi Agarwal

Description:
    For a vertex labelled graph, output a heat map pdf file which shows
    the frequency of connections between any two labels in the graph
    Plot is saved at the supplied pdf path titled with the supplied plot title

USAGE: 
    python label_connection_confusion_matrix.py vertex_label_file edge_file output_pdf_path plot_title
"""

###############################################################################

import sys
import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import seaborn as sns
import math

if len(sys.argv) != 5:
    print("""
Description:
    For a vertex labelled graph, output a heat map pdf file which shows
    the frequency of connections between any two labels in the graph
    Plot is saved at the supplied pdf path titled with the supplied plot title

USAGE: 
    python label_connection_confusion_matrix.py vertex_label_file edge_file output_pdf_path plot_title
"""
)
    
def format_with_k_m(value, _):
    if value >= 1e3:
        return f'{value / 1e3:.0f}K'
    else:
        return str(value)

# Read vertex labels
with open(sys.argv[1]) as f:
    raw = f.readlines()

vertices = dict()
labels = set()

for x in raw:
    v, l = x.strip().split()
    vertices[v] = l
    labels.add(l)

labels = list(labels)

# print(f"{len(vertices)} {len(labels)}")

# Read edges
with open(sys.argv[2]) as f:
    raw = f.readlines()

l2l_edges = np.zeros((len(labels), len(labels)))

counter = 0
for e in raw:
    u, v = e.strip().split()[:2]
    pos_u = labels.index(vertices[u])
    pos_v = labels.index(vertices[v])
    l2l_edges[pos_u, pos_v] = l2l_edges[pos_u, pos_v] + 1
    l2l_edges[pos_v, pos_u] = l2l_edges[pos_v, pos_u] + 1
    counter = counter + 1
    if counter%1000000==0:
        print(f"{counter}", end=' ', flush=True)
print("done\nPreparing heatmap")

# print(f"{l2l_edges}")

tick_size = math.floor(len(labels)/8)

df = DataFrame(l2l_edges, index=range(1, len(labels)+1), columns=range(1, len(labels)+1))
# print(df)

sns.set(font_scale=1.65)
s=sns.heatmap(df, cmap='rocket_r', yticklabels=tick_size, xticklabels=tick_size)

# # Set the font size for the x-axis and y-axis labels
# s.set_xticklabels(s.get_xticklabels(), fontsize=12)
# s.set_yticklabels(s.get_yticklabels(), fontsize=12)

# cbar = s.collections[0].colorbar
# cbar.formatter = FuncFormatter(format_with_k_m)
# cbar.update_ticks()

# plt.title(sys.argv[4], fontsize=14, fontweight="bold")
plt.savefig(sys.argv[3], format="pdf", bbox_inches="tight")