import os
import pandas as pd
import numpy as np

# Grab all the PubMed IDs for a particular keyword and check headlines/abstracts for gene names
# Must have EDirect installed on command line and api key specified in .zshrc

term = input("Search PubMed for term: ")
term = "+".join(term.split(" "))
term = term.lower()

# Get articles
for year in range(2000,2026):
	print(f"Searching \"{term}\" from {year}...")
	data = os.popen(f"esearch -db pubmed -query \"{term}\" -mindate {year} -maxdate {year} | efetch -format abstract").read()
	with open(f"data/articles/pubmed_{term}_{year}.txt", "w") as f:
		f.write(data)

