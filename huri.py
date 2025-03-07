import pandas as pd

# HuRI is a list of protein-protein interactions

def huri(raw=False):
	if raw:
		return pd.read_csv("data/HuRI.tsv", sep = "\t", header = None)
	else:
		return pd.read_csv("data/HuRI_genenames.tsv", sep = "\t", header = None)



