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

# Pass in pandas DF of pairwise distances:
# path_pairs = pd.DataFrame(nx.all_pairs_shortest_path_length(g))
def avg_path_length(path_pairs):
	dists = np.zeros((len(path_pairs),len(path_pairs)))
	for i in range(len(path_pairs)):
		dists[i,:] = list(path_pairs.iloc[i,1].values())
	
	dists = dists.flatten()
	dists = dists[dists != 0]
	
	return dists.mean(), dists.std()
	#return f"{round(dists.mean(),2)}±{round(dists.std(),2)}"

# Pass in pairwise distances and Series or list of genes of interest
def avg_path_length_genes(path_pairs, genes_of_interest):
	genes_of_interest = pd.Series([x for x in genes_of_interest if x in path_pairs.iloc[0,1].keys()]) # Only include genes present in the graph
	path_pairs = path_pairs[path_pairs.iloc[:,0].isin(genes_of_interest)]
	
	dists = np.zeros((len(path_pairs),len(path_pairs)))
	for i in tqdm(range(len(path_pairs))):
		dists[i,:] = pd.DataFrame(path_pairs.iloc[i,1],index=range(1)).loc[0,genes_of_interest]
	
	dists = dists.flatten()
	dists = dists[dists != 0]
	
	return dists.mean(), dists.std()
	#return f"{round(dists.mean(),2)}±{round(dists.std(),2)}"


pairs = pd.read_csv("data/HuRI_genenames.tsv", sep = "\t", header = None)
genes = list(set(list(pairs.iloc[:,0])+list(pairs.iloc[:,1])))
genes.sort()

g = nx.Graph()
g.add_nodes_from(genes)
for i in range(len(pairs)):
	g.add_edge(pairs.iloc[i,0], pairs.iloc[i,1])

# There are 73 disconnected subgraphs. All but the first one have at most 3 nodes. First graph contains 8109 nodes.
# g = [g.subgraph(c).copy() for c in nx.connected_components(g)]
# for i in range(4,8):
# 	plt.subplot(2,2,i-3)
# 	layout = nx.spring_layout(g[i], k=0.1, scale=0.01, iterations=1)
# 	nx.draw(g[i], layout, with_labels=True, node_size=100, node_color='skyblue', font_weight='bold')
# plt.show()

g = remove_disconnected(g)
path_pairs = pd.DataFrame(nx.all_pairs_shortest_path_length(g))
#print(avg_path_length(path_pairs)) # 3.84±0.91

files = os.listdir("data/results")
files.sort()
disorder = []
path_length_avgs = []
path_length_stds = []
for f in files:
	print(f)
	disorder.append(" ".join(f.split("_")[1].split("+")))
	genes_of_interest = pd.read_csv(f"data/results/{f}", sep = "\t", header = None)
	#genes_of_interest = genes_of_interest[genes_of_interest.iloc[:,1] > 1]
	path_length_avg, path_length_std = avg_path_length_genes(path_pairs,genes_of_interest.iloc[:,0])
	path_length_avgs.append(path_length_avg)
	path_length_stds.append(path_length_std)

plt.bar(disorder, path_length_avgs, yerr=path_length_stds, capsize=5, edgecolor="black")
plt.xticks(rotation=90, fontsize=8)
plt.axhline(y=3.84, color='black', linestyle='--', linewidth=2)
plt.ylabel("Average path length between associated genes")
plt.tight_layout()
plt.show()








