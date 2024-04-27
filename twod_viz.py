import sys
import networkx as nx
import matplotlib.pyplot as plt

colordict = {"A": "lightblue", "U": "red", "G": "lightgreen", "C": "yellow"}
G = nx.Graph()
rna_seq = sys.argv[1]
# rna_seq = "GAUGCUGGUC"
node_labels = []
color_seq = []
for i in range(len(rna_seq)):
    node_labels.append(rna_seq[i] + str(i))
    color_seq.append(colordict[rna_seq[i]])
edges_list = []
for i in range(len(rna_seq) - 1):
    edges_list.append((rna_seq[i] + str(i), rna_seq[i + 1] + str(i + 1)))
# print(node_labels)
G.add_nodes_from(node_labels)
print(edges_list)
csvstr = sys.argv[2]
str_arr = csvstr.split(",")
for i in range(int(len(str_arr) / 2)):
    edges_list.append(
        (
            rna_seq[int(str_arr[2 * i])] + str_arr[2 * i],
            rna_seq[int(str_arr[2 * i + 1])] + str_arr[2 * i + 1],
        )
    )
G.add_edges_from(edges_list)

pos = nx.circular_layout(G)
nx.draw_networkx(G, with_labels=True, node_color=color_seq)
plt.axis("off")
plt.show()
