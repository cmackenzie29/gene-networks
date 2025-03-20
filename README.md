# Network Analysis of the Human Reference Interactome
[Read the full report](report.pdf)
## Results
Preliminary analysis of the HuRI network revealed 73 disconnected subgraphs, with one of these (the
“main network”) containing 8,109 genes and each of the 72 others containing no more than three genes
(Fig. 1A-B). For simplicity, I discarded the small subgraphs and considered only the features of the main
network. Degree-rank plots revealed that while more than half the genes in the network have fewer than
10 connections, some have several hundred connections (Fig. 1C-D). I also analyzed the betweenness
centrality of each node, which was correlated with the node’s degree with a Pearson’s coefficient of 0.76
(Fig. 1E).

---

![Image](/figures/Fig1.png)

Figure 1. (A) Visualization of HuRI protein interaction network, with node diameter representing degree
and color representing betweenness centrality (red: higher, blue: lower). (B) Examples of small,
disconnected sub-graphs present in the HuRI dataset. (C) Degree-rank plot of all nodes in the main
network, with degree plotted on a logarithmic scale. (D) Plot showing the number of nodes of each degree
in the main network, plotted on a log-log scale. (E) Scatterplot showing the relationship between each
node’s degree and betweenness centrality, plotted on a log-log scale.

---

Next, I analyzed autism-specific genes within the network. I used PubMed to identify 536 autism-related genes within the main network, spanning a wide range of centralities (Fig. 2-left). While I had
suspected that the most-referenced genes would likely correspond to those with more central roles in the
network, I found essentially no correlation between degree and number of PubMed abstracts referencing
the gene (Fig. 2-right). I also expected to find that the autism-related genes might be clustered together in
the network, so I was surprised to find that the mean ± S.D. path length between autism-related genes was
3.89±0.93, marginally higher than the average path length between all genes in the network, 3.84±0.91. I
found this to be the case for 34 other disorders as well, with the average path length between disease-
related genes differing very minimally from the average path length between all genes (Supp. Fig. 1). A
summary of the PubMed search results for all 35 disorders considered is given in Table 1.

---

![Image](/figures/Fig2.png)

Figure 2. (Left) Visualization of the network with autism-related genes emphasized. (Right) Scatterplot
showing the relationship between a gene’s degree in the network and the number of times the gene was
referenced in PubMed abstracts since 2000 matching the search term “autism.”

![Image](/figures/FigS1.png)

Supplemental Figure 1. Average path length between each set of disease-relevant genes. Each bar
represents the mean ± standard deviation.

---

Table 1. Summary of PubMed search results for each disorder

| Disorder | Article Count | Genes in Network |
| ---------- | ---------- | ---------- |
| Alzheimer's | 234,989 | 1,075 |
| Anxiety | 344,106 | 521 |
| Arthritis | 308,327 | 1,135 |
| ADHD | 54,383 | 187 |
| Autism | 86,838 | 536 |
| Bipolar | 81,368 | 445 |
| Celiac | 27,962 | 137 |
| Chronic fatigue syndrome | 9,626 | 67 |
| Congenital heart disease | 137,426 | 508 |
| Crohn's disease | 58,671 | 387 | 387 |
| Depression | 558,886 | 902 |
| Diabetes | 906,510 | 2,132 |
| Epilepsy | 144,985 | 796 |
| Fibromyalgia | 14,035 | 61 |
| Hypercholesterolemia | 36,668 | 280 |
| Hyperthyroidism | 24,621 | 190 |
| Hypothyroidism | 32,787 | 233 |
| Idiopathic pulmonary fibrosis | 16,885 | 276 |
| Lupus erythematosus | 60,294 | 468 |
| Melanoma | 135,069 | 1,220 |
| Meniere's | 4,884 | 36 |
| Migraine | 40,428 | 153 |
| Multiple sclerosis | 103,020 | 567 |
| Myasthenia gravis | 12,387 | 103 |
| Parkinson's | 166,341 | 773 |
| PCOS | 19,064 | 260 |
| Peripheral neuropathy | 138,204 | 716 |
| Postural orthostatic tachycardia | 1,460 | 11 |
| Psoriasis | 54,692 | 440 |
| Schizophrenia | 130,431 | 615 |
| Sjogren's | 18,672 | 224 |
| Substance use disorder | 228,189 | 439 |
| Tourette syndrome | 5,363 | 40 |
| Transverse myelitis | 7,677 | 53 |
| Ulcerative colitis | 49,963 | 413 |
| **Total** | 4,255,211 | 8,109 |

---

Next, I wanted to understand the extent to which the removal of one or a few nodes might disrupt
the broader connectivity of the network. I started by removing the highest-degree node, CYSRT1, with
498 connections, which slightly increased the average path length of the network to 3.86±0.91. Next, I
found that removing the most-referenced autism-related gene, FMR1, which only had 9 connections, also
increased the average path length between autism-related genes to 3.90±0.93. While these increases were
very subtle, I was surprised that such small changes in the network could have a detectable effect at all.

I extended this perturbation analysis by examining the effects of removing the single most-
referenced gene and the twenty most-referenced genes (together) for each disorder. While these
perturbations had almost no effect on the average path length between all genes in the network, they did
tend to increase the average path length between each set of corresponding disease-related genes (Fig. 3).
Specifically, removing the single most-referenced gene increased the average disease-related gene path
length by 0.00015±0.00508 (not significant, p>0.05, by t-test comparison to 0 mean), while removal of
the twenty most-referenced gene increased the average disease-related gene path length by
0.02137±0.053046 (significant, p<0.05, by t-test comparison to 0 mean). In short, these results suggest
that perturbation of a disease-specific gene may affect the connectivity of the genes associated with the
same disease more than connectivity of the genes in the broader network.

---

![Image](/figures/Fig3.png)

Figure 3. (Top) Comparison of average path lengths between disease-related genes (blue) vs. path lengths
between these genes following the removal of the single most-referenced (orange) or twenty most-
referenced (green) genes for each disorder. (Bottom) Comparison of average path lengths among all genes
in the network following the removal of the single most-referenced (blue) or twenty most-referenced
(orange) genes for each disorder.

## Resources
- **Gene names**: [https://www.genenames.org/download/custom/](https://www.genenames.org/download/custom/)
- **Protein-Protein Interactions (HuRI)**: Luck K, Kim DK, Lambourne L, et al. A reference map of the human binary protein interactome. Nature. 2020;580(7803):402-408. [doi:10.1038/s41586-020-2188-x](https://pmc.ncbi.nlm.nih.gov/articles/PMC7169983/)

