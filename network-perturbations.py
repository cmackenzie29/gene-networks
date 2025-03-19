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

def avg_path_length(path_pairs):
	dists = np.zeros((len(path_pairs),len(path_pairs)))
	for i in range(len(path_pairs)):
		dists[i,:] = list(path_pairs.iloc[i,1].values())
	
	dists = dists.flatten()
	dists = dists[dists != 0]
	
	return dists.mean(), dists.std()

def avg_path_length_genes(path_pairs, genes_of_interest):
	genes_of_interest = pd.Series([x for x in genes_of_interest if x in path_pairs.iloc[0,1].keys()]) # Only include genes present in the graph
	path_pairs = path_pairs[path_pairs.iloc[:,0].isin(genes_of_interest)]
	
	dists = np.zeros((len(path_pairs),len(path_pairs)))
	for i in tqdm(range(len(path_pairs))):
		dists[i,:] = pd.DataFrame(path_pairs.iloc[i,1],index=range(1)).loc[0,genes_of_interest]
	
	dists = dists.flatten()
	dists = dists[dists != 0]
	
	return dists.mean(), dists.std()

pairs = pd.read_csv("data/HuRI_genenames.tsv", sep = "\t", header = None)
genes = list(set(list(pairs.iloc[:,0])+list(pairs.iloc[:,1])))
genes.sort()

g = nx.Graph()
g.add_nodes_from(genes)
for i in range(len(pairs)):
	g.add_edge(pairs.iloc[i,0], pairs.iloc[i,1])

g = remove_disconnected(g)

# degree_sequence = sorted((d for n, d in g.degree()), reverse=True)
# dmax = max(degree_sequence)

# Explore effect of dropping CYSRT1, which has the highest degree in the network (498)
g_drop = g.copy()
g_drop.remove_node("CYSRT1")
g_drop = remove_disconnected(g_drop)
print(f"{len(g.nodes) - len(g_drop.nodes) - 1} additional genes excluded from network") # 10
path_pairs = pd.DataFrame(nx.all_pairs_shortest_path_length(g_drop))
print(avg_path_length(path_pairs)) # 3.86±0.91

# Explore effect of dropping top autism-related gene (FMR1)
autism_genes = pd.read_csv(f"data/results/pubmed_autism_results.txt", sep = "\t", header = None)
autism_genes = autism_genes.sort_values(by=1)

g_drop = g.copy()
g_drop.remove_node("FMR1") # Interestingly, FMR1 is only degree 9
g_drop = remove_disconnected(g_drop)
print(f"{len(g.nodes) - len(g_drop.nodes) - 1} additional genes excluded from network") # 1
path_pairs = pd.DataFrame(nx.all_pairs_shortest_path_length(g_drop))
print(avg_path_length(path_pairs)) # 3.84±0.91
print(avg_path_length_genes(path_pairs, autism_genes.iloc[:,0])) # 3.90±0.93

# Explore effect of dropping top 20 autism-related genes
g_drop = g.copy()
for gene in autism_genes.iloc[-20:,0]:
	g_drop.remove_node(gene)

g_drop = remove_disconnected(g_drop)
print(f"{len(g.nodes) - len(g_drop.nodes) - 19} additional genes excluded from network") # 10
path_pairs = pd.DataFrame(nx.all_pairs_shortest_path_length(g_drop))
print(avg_path_length(path_pairs)) # 3.84±0.91
print(avg_path_length_genes(path_pairs, autism_genes.iloc[:,0])) # 3.91±0.93

















