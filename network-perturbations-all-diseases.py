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

# Explore effect of dropping top gene and top 20 genes for each disease
disorder = []
none_dropped_path_length_genes = []
one_dropped_genes_excluded = []
one_dropped_path_length_all = []
one_dropped_path_length_genes = []
twenty_dropped_genes_excluded = []
twenty_dropped_path_length_all = []
twenty_dropped_path_length_genes = []

complete_network_path_pairs = pd.DataFrame(nx.all_pairs_shortest_path_length(g))

files = os.listdir("data/results")
files.sort()
for f in files:
	print(f)
	genes_of_interest = pd.read_csv(f"data/results/{f}", sep = "\t", header = None)
	genes_of_interest = genes_of_interest[genes_of_interest.iloc[:,0].isin(g.nodes())] # Remove genes not in the network
	genes_of_interest = genes_of_interest.sort_values(by=1)

	if len(genes_of_interest) < 25:
		continue

	disorder.append(" ".join(f.split("_")[1].split("+")))
	none_dropped_path_length_genes.append(avg_path_length_genes(complete_network_path_pairs, genes_of_interest.iloc[:,0]))

	# Effect of dropping top gene
	g_drop = g.copy()
	g_drop.remove_node(genes_of_interest.iloc[-1,0])
	g_drop = remove_disconnected(g_drop)
	one_dropped_genes_excluded.append(len(g.nodes) - len(g_drop.nodes) - 1)
	path_pairs = pd.DataFrame(nx.all_pairs_shortest_path_length(g_drop))
	one_dropped_path_length_all.append(avg_path_length(path_pairs))
	one_dropped_path_length_genes.append(avg_path_length_genes(path_pairs, genes_of_interest.iloc[:,0]))

	# Effect of dropping top 20 genes
	g_drop = g.copy()
	for gene in genes_of_interest.iloc[-20:,0]:
		g_drop.remove_node(gene)

	g_drop = remove_disconnected(g_drop)
	twenty_dropped_genes_excluded.append(len(g.nodes) - len(g_drop.nodes) - 20)
	path_pairs = pd.DataFrame(nx.all_pairs_shortest_path_length(g_drop))
	twenty_dropped_path_length_all.append(avg_path_length(path_pairs))
	twenty_dropped_path_length_genes.append(avg_path_length_genes(path_pairs, genes_of_interest.iloc[:,0]))

# Plot number of genes dropped from the network
n_bars = len(disorder)
bar_width = 0.35
index = np.arange(n_bars)
plt.bar(index, one_dropped_genes_excluded, bar_width, label='Top gene removed')
plt.bar(index + bar_width, twenty_dropped_genes_excluded, bar_width, label='Top 20 genes removed')
plt.ylabel('# of additional genes removed')
plt.xticks(index + bar_width / 2, disorder, rotation=90, fontsize=8)
plt.legend()
plt.tight_layout()
plt.show()

# Plot average path lengths among disorder genes in the network
n_bars = len(disorder)
bar_width = 0.35
index = np.arange(n_bars) * 1.5
plt.bar(index, [x[0] for x in none_dropped_path_length_genes], bar_width, label='Complete network')
plt.bar(index + bar_width, [x[0] for x in one_dropped_path_length_genes], bar_width, label='Top gene removed')
plt.bar(index + 2*bar_width, [x[0] for x in twenty_dropped_path_length_genes], bar_width, label='Top 20 genes removed')
plt.ylabel('# of additional genes removed')
plt.xticks(index + bar_width, disorder, rotation=90, fontsize=8)
plt.legend()
plt.tight_layout()
plt.show()

# Plot average path lengths in the whole network
n_bars = len(disorder)
bar_width = 0.35
index = np.arange(n_bars)
plt.bar(index, [x[0] for x in one_dropped_path_length_all], bar_width, label='Top gene removed')
plt.bar(index + bar_width, [x[0] for x in twenty_dropped_path_length_all], bar_width, label='Top 20 genes removed')
plt.axhline(y=3.84397798, color='black', linestyle='--', linewidth=1.0)
plt.ylabel('# of additional genes removed')
plt.xticks(index + bar_width / 2, disorder, rotation=90, fontsize=8)
plt.legend()
plt.tight_layout()
plt.show()


# Cache results
# import pickle
# with open("/data/gene_perturbation_cache/disorder.sav", "wb") as f:
# 	pickle.dump(disorder, f)

# with open("/data/gene_perturbation_cache/none_dropped_path_length_genes.sav", "wb") as f:
# 	pickle.dump(none_dropped_path_length_genes, f)

# with open("/data/gene_perturbation_cache/one_dropped_genes_excluded.sav", "wb") as f:
# 	pickle.dump(one_dropped_genes_excluded, f)

# with open("/data/gene_perturbation_cache/one_dropped_path_length_all.sav", "wb") as f:
# 	pickle.dump(one_dropped_path_length_all, f)

# with open("/data/gene_perturbation_cache/one_dropped_path_length_genes.sav", "wb") as f:
# 	pickle.dump(one_dropped_path_length_genes, f)

# with open("/data/gene_perturbation_cache/twenty_dropped_genes_excluded.sav", "wb") as f:
# 	pickle.dump(twenty_dropped_genes_excluded, f)

# with open("/data/gene_perturbation_cache/twenty_dropped_path_length_all.sav", "wb") as f:
# 	pickle.dump(twenty_dropped_path_length_all, f)

# with open("/data/gene_perturbation_cache/twenty_dropped_path_length_genes.sav", "wb") as f:
# 	pickle.dump(twenty_dropped_path_length_genes, f)











