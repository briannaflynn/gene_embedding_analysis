#!/usr/bin/env python3
#Date: Sept 3, 2024; Oct 4, 2024; Oct 14, 2024
#Author: Muyoung Lee
#Description:
#1. From the reference human proteome genes, gather genes belonging to the input feature table
#2. Calculate Spearman correltaions between all gene pairs.
#3. Sort gene pairs based on calculated correlations.
#4. If the two genes in a gene pair share at least one GOBP terms whose level is larger than 2 (deeper than 2), mark this gene pair as positive. Consider ancestors GOBP terms as well.
#5. Based on step 3 and 4, calculate AUROC and AUPRC.

import pickle
import numpy as np
import pandas as pd
import random
from sklearn import metrics

def cal_AUROC(sorted_gene_list, positive_set):
	total = len(sorted_gene_list)
	actual_positive = positive_set & set(sorted_gene_list)
	dy_TPR = 1 / len(actual_positive)
	dx_FPR = 1 / (total - len(actual_positive))
	x, y, AUROC = 0, 0, 0
	for gene in sorted_gene_list:
		dAUROC = 0
		if gene in actual_positive:
			y += 1
		else:
			x += 1
			dAUROC = dx_FPR * y * dy_TPR
		AUROC += dAUROC
	remainder_AUROC = (1 - x * dx_FPR) * (y * dy_TPR + 1) * 0.5
	total_AUROC = AUROC + remainder_AUROC
	return total_AUROC

def cal_AUPRC(sorted_gene_list, positive_set):
	total = len(sorted_gene_list)
	actual_positive = positive_set & set(sorted_gene_list)
	# x: recall = (TP / Actual Positive)
	# y: precision = (TP / Predicted Positive)
	AP = len(actual_positive)
	PP = 0
	TP = 0
	calculated_values = []
	for gene in sorted_gene_list:
		PP += 1
		if gene in actual_positive:
			TP += 1
		calculated_values.append((TP/AP, TP/PP))
	recall = [x[0] for x in calculated_values]
	precision = [x[1] for x in calculated_values]
	AUC = metrics.auc(recall, precision)
	return AUC

# GOBP_level: GOBP -> level
with open("GOBP/GOBP_level.pkl", "rb") as INPUT:
	GOBP_level = pickle.load(INPUT)

# GOBP_BPancestors: GOBP -> Ancestor GOBP terms
with open("GOBP/GOBP_BP_ancestors.pkl", "rb") as INPUT2:
	GOBP_BPancestors = pickle.load(INPUT2)

ensID_BPterms = {}
with open("GOBP/EnsID_GOid_GOterm_GOdomain.BP_only.tsv") as INPUT:
	for line in INPUT:
		words = line.strip().split("\t")
		ensID, goid = words[0], words[1]
		if goid not in GOBP_level:
			continue
		elif GOBP_level[goid] <= 2:
			continue			
		if ensID not in ensID_BPterms:
			ensID_BPterms[ensID] = set([goid])
		else:
			ensID_BPterms[ensID].add(goid)
		if goid in GOBP_BPancestors:
			ancestors = []
			for GOBP in GOBP_BPancestors[goid]:
				if GOBP_level[GOBP] > 2:
					ancestors.append(GOBP)
			if len(ancestors) > 0:
				ensID_BPterms[ensID].update(ancestors)

reference_proteome = set()
with open("human_reference_proteome/HumanRefProteome_UniProtACC_EnsGeneID.tsv") as INPUT:
	for line in INPUT:
		ensID = line.strip().split("\t")[1]
		if ensID != "N/A":
			ensID_list = ensID.split(",")
			reference_proteome.update(ensID_list)

# output from CZ_CellxGene/make_a_tidy_table.py
df = pd.read_csv("CZ_CellxGene/ratio_of_n_and_n_cells_cell_type.gene.tissue-celltype.tsv", sep="\t", header=0, index_col=0)

domain = list(set(df.index) & reference_proteome)
print(len(domain))

#sampled_genes = random.sample(domain, 1000)
# Now, we use all genes.
sampled_genes = domain

gene_gene2_corr = {}
na_list = []
for gene in sampled_genes:
	for gene2 in domain:
		if gene2 == gene:
			continue
		gene, gene2 = sorted([gene, gene2])
		if (gene, gene2) in gene_gene2_corr:
			continue
		elif (gene, gene2) in na_list:
			continue
		else:
			cor_value = df.loc[gene].corr(df.loc[gene2], method="spearman", min_periods=30)
			if np.isnan(cor_value):
				na_list.append((gene, gene2))
			else:
				gene_gene2_corr[(gene,gene2)] = cor_value

gene_sorted_corr = dict(sorted(gene_gene2_corr.items(), key=lambda item: item[1], reverse=True))
del gene_gene2_corr

for (gene1, gene2) in gene_sorted_corr:
	positive_set = set() # "gene2" that shares at least 1 GOBP (level>2) with "gene"
	if gene1 in ensID_BPterms:
		gene1_BP_set = ensID_BPterms[gene]
	else:
		gene1_BP_set = set()
	
	if gene2 in ensID_BPterms:
		gene2_BP_set = ensID_BPterms[gene2]
	else:
		gene2_BP_set = set()

	count = len(gene1_BP_set & gene2_BP_set)
	if count > 0:
		positive_set.add((gene1, gene2))

	if len(positive_set) > 0:
		AUROC = cal_AUROC(list(gene_sorted_corr) + na_list, positive_set)
		AUPRC = cal_AUPRC(list(gene_sorted_corr) + na_list, positive_set)
	else:
		AUROC = np.nan
		AUPRC = np.nan

print(f"AUROC: {AUROC}, AUPRC: {AUPRC}")