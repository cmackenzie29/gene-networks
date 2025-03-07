import pandas as pd

# HuRI is a list of protein-protein interactions
# Script converts Ensembl gene IDs to their names

ensembl = pd.read_csv("data/ensembl.tsv", sep = "\t")
huri = pd.read_csv("data/HuRI.tsv", sep = "\t", header = None)

for i in range(len(huri)-1, -1, -1):
	print(i)
	for j in range(2):
		gene_name = ensembl[ensembl.iloc[:,2]==huri.iloc[i,j]]
		if len(gene_name) > 0:
			huri.loc[i,j] = gene_name.iloc[0,0]
		else:
			huri = huri.drop(i)
			break

huri.to_csv("data/HuRI_genenames.tsv", index = False, header = False, sep = "\t")

