import os
import sys
from tqdm import tqdm
import pandas as pd
import numpy as np

# Grab all the PubMed IDs for a particular keyword and check headlines/abstracts for gene names
# Must have EDirect installed on command line and api key specified in .zshrc

term = input("Search PubMed for term: ")
term = "+".join(term.split(" "))
term = term.lower()

if os.path.isfile(f"data/{term}_pubmed_results.txt"): # Results already exist for term and will be overwritten.
	overwrite_permission = input(f"\"{term}\" has saved results that appear to be final. Overwrite? (y/n): ")
	if overwrite_permission != "y":
		sys.exit()

print(f"Searching \"{term}\"...")

gene_names = pd.read_csv("data/ensembl.tsv", sep = "\t")
gene_names = gene_names.apply(lambda x: x.str.lower())
gene_names = dict(zip(gene_names.iloc[:,0], gene_names.iloc[:,1]))
genes = pd.read_csv("data/HuRI_genenames.tsv", sep = "\t")
genes = list(set(list(genes.iloc[:,0])+list(genes.iloc[:,1])))
for g in range(len(genes)):
	genes[g] = genes[g].lower()
genes.sort()
hits = np.zeros(len(genes)).astype(int)
num_articles = 0

# Get articles
for year in range(2000,2026):
	data = os.popen(f"esearch -db pubmed -query \"{term}\" -mindate {year} -maxdate {year} | efetch -format abstract").read()
	data = data.split("\n\n\n")
	data = np.array(data)
	data = np.char.replace(data, " \n", " ")
	data = np.char.replace(data, "\n ", " ")
	data = np.char.replace(data, "\n", "")
	data = np.char.lower(data)

	print(f"Analyzing {len(data)} articles from {year}")
	num_articles += len(data)

	for g in tqdm(range(len(genes))):
		hits[g] += np.logical_or.reduce(((np.char.find(data, f"{genes[g]} expression") >= 0), (np.char.find(data, f"{genes[g]} gene") >= 0), (np.char.find(data, f"gene {genes[g]}") >= 0), (np.char.find(data, gene_names[genes[g]]) >= 0))).sum()

print(f"{num_articles} articles analyzed")
results = pd.DataFrame({"A": genes, "B": hits})
results = results[results.B>0]
if len(results) > 0:
	results.A = results.A.apply(lambda x: x.upper())
	results.to_csv(f"data/{term}_results.tsv", index = False, header = False, sep = "\t")
else:
	print("Analysis finished with no matches")

























