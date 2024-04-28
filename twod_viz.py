"""
@file twod_viz.py
@brief This file contains a Python script that generates a two-dimensional visualization of the secondary structure of an RNA sequence using NetworkX and Matplotlib.

The script takes two command-line arguments:
1. The RNA sequence as a string.
2. A comma-separated string representing the bonded pairs of indices in the RNA sequence.

The script creates a graph object using NetworkX, with nodes representing the nucleotides in the RNA sequence and edges representing the bonds between nucleotides in the secondary structure. The nodes are colored according to the nucleotide type (A, U, G, or C). The graph is then visualized using Matplotlib's circular layout.

@author Y.S Sandeep
@author Suman Sekhar Sahoo
@author Raghuram Venkatesan
@author Bhaarat K
@author Kamalesh Ram
"""

import sys
import networkx as nx
import matplotlib.pyplot as plt

"""
@var global colordict
A dictionary that maps nucleotide types to colors for visualization.
"""
colordict = {"A": "lightblue", "U": "red", "G": "lightgreen", "C": "yellow"}

"""
@var global G
The NetworkX graph object.
"""
G = nx.Graph()

"""
@var global rna_seq
The RNA sequence as a string.
"""
rna_seq = sys.argv[1]

"""
@var global node_labels
A list of node labels, where each label is a nucleotide type concatenated with its index in the sequence.
"""
node_labels = []

"""
@var global color_seq
A list of node colors corresponding to the nucleotide types.
"""
color_seq = []
edge_color_seq = []
line_widths = []
"""
@var global edges_list
A list of edges representing the bonds between adjacent nucleotides in the sequence.
"""
edges_list = []

for i in range(len(rna_seq)):
    node_labels.append(rna_seq[i] + str(i))
    color_seq.append(colordict[rna_seq[i]])
for i in range(len(rna_seq) - 1):
    edges_list.append((rna_seq[i] + str(i), rna_seq[i + 1] + str(i + 1)))
    edge_color_seq.append("black")
    line_widths.append(1.0)
G.add_nodes_from(node_labels)
# print(edges_list)

"""
@var global csvstr
The comma-separated string representing the bonded pairs of indices.
"""
csvstr = sys.argv[2]

"""
@var global str_arr
An array of indices obtained by splitting the `csvstr` string.
"""
str_arr = csvstr.split(",")

for i in range(int(len(str_arr) / 2)):
    edges_list.append(
        (
            rna_seq[int(str_arr[2 * i])] + str_arr[2 * i],
            rna_seq[int(str_arr[2 * i + 1])] + str_arr[2 * i + 1],
        )
    )
    edge_color_seq.append("red")
    line_widths.append(1.2)
G.add_edges_from(edges_list)

"""
@var global pos
A dictionary mapping nodes to their positions in the visualization, obtained using NetworkX's circular layout.
"""
pos = nx.kamada_kawai_layout(G)

"""
@var global with_labels
A boolean flag indicating whether to display node labels in the visualization.
"""

"""
@var global True
A boolean value set to True, indicating that node labels should be displayed.
"""

"""
@var global node_color
A list of colors for the nodes, corresponding to the nucleotide types.
"""
# nx.draw_networkx(G, pos, with_labels=True, node_color=color_seq)
nx.draw_networkx_nodes(G, pos, node_size=300, node_color=color_seq)
nx.draw_networkx_edges(
    G, pos, edgelist=edges_list, width=line_widths, edge_color=edge_color_seq
)
nx.draw_networkx_labels(G, pos, font_size=10)
plt.axis("off")
plt.show()
