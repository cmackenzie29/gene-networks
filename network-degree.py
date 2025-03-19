import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as cols
import pandas as pd
import numpy as np
import os
from tqdm import tqdm

def remove_disconnected(g):
	sub_g = [g.subgraph(c).copy() for c in nx.connected_components(g)]
	sub_sizes = np.array([len(x.nodes) for x in sub_g])
	return sub_g[sub_sizes.argmax()]

pairs = pd.read_csv("data/HuRI_genenames.tsv", sep = "\t", header = None)
genes = list(set(list(pairs.iloc[:,0])+list(pairs.iloc[:,1])))
genes.sort()

g = nx.Graph()
g.add_nodes_from(genes)
for i in range(len(pairs)):
	g.add_edge(pairs.iloc[i,0], pairs.iloc[i,1])

g = remove_disconnected(g)

degree_sequence = sorted((d for n, d in g.degree()), reverse=True)
dmax = max(degree_sequence)

plt.subplot(1,2,1)
plt.plot(np.log(degree_sequence), marker="o")
plt.title("A", loc="left", fontsize=20)
plt.ylabel("log(Degree)", fontsize=15)
plt.xlabel("Rank", fontsize=15)

plt.subplot(1,2,2)
plt.scatter(*np.log(np.unique(degree_sequence, return_counts=True)))
plt.title("B", loc="left", fontsize=20)
plt.xlabel("log(Degree)", fontsize=15)
plt.ylabel("log(# of Nodes)", fontsize=15)
plt.show()


# Best: k = 0.15
layout = nx.spring_layout(g, k=0.15, seed=1)
centrality = nx.betweenness_centrality(g)
max_centrality = max(list(centrality.values()))

color_map = []
for node in g.nodes:
	c = centrality[node]/max_centrality
	color_map.append(plt.cm.jet((101*c)/(1+100*c)))

node_sizes = [g.degree(node)**0.8 for node in g.nodes()]
plt.figure(figsize=(7, 7))
nx.draw_networkx_nodes(g, layout, node_size=node_sizes, node_color=color_map, edgecolors='black', linewidths=0.1)
nx.draw_networkx_edges(g, layout, node_size=node_sizes, width=0.2, edge_color=cols.colorConverter.to_rgba('black', alpha=0.1))
plt.show()

degree = []
betweenness = []
for node in g.degree():
	degree.append(node[1])
	betweenness.append(centrality[node[0]])
plt.scatter(np.log(degree), np.log(betweenness), s=3.0)
plt.xlabel("log(Degree)", fontsize=15)
plt.ylabel("log(Betweenness Centrality)", fontsize=15)
plt.show()

# np.corrcoef(degree, betweenness) # Pearson correlation: 0.76



# Highlight the autism-related genes in the network
genes_of_interest = pd.read_csv(f"data/results/pubmed_autism_results.txt", sep = "\t", header = None).iloc[:,0]
alpha_map = np.zeros(0)
for node in g.nodes:
	if (genes_of_interest == node).sum() > 0:
		alpha_map = np.append(alpha_map, 1.0)
	else:
		alpha_map = np.append(alpha_map, 0.05)

plt.figure(figsize=(7, 7))
nx.draw_networkx_nodes(g, layout, node_size=node_sizes, node_color=color_map, edgecolors='black', linewidths=0.1, alpha=alpha_map)
nx.draw_networkx_edges(g, layout, node_size=node_sizes, width=0.2, edge_color=cols.colorConverter.to_rgba('black', alpha=0.1))
plt.show()


# Examine relationship between most cited autism genes and their degree in the network
genes_of_interest = pd.read_csv(f"data/results/pubmed_autism_results.txt", sep = "\t", header = None)
degree = []
pm_references = []
for node in g.degree():
	if (genes_of_interest.iloc[:,0] == node[0]).sum() > 0:
		degree.append(node[1])
		pm_references.append(genes_of_interest[genes_of_interest.iloc[:,0] == node[0]].iloc[0,1])

plt.scatter(np.log(degree), np.log(pm_references), s=3.0)
plt.xlabel("log(Degree)", fontsize=15)
plt.ylabel("log(# Pubmed References)", fontsize=15)
plt.show()





