#!/usr/bin/env python3
#Date: Sept 3, 2024; Oct 4, 2024; Oct 14, 2024; Oct 23, 2024
#Author: Muyoung Lee
#Description:
#1. From the reference human proteome genes, gather genes belonging to the input feature table.
#2. Calculate Spearman correlations between all gene pairs.
#3. Sort gene pairs based on calculated correlations.
#4. If the two genes in one gene pair share at least one GOBP term whose level is larger than 2 (deeper than 2), mark this gene pair as positive. Consider ancestor GOBP terms as well.
#5. Based on steps 3 and 4, calculate AUROC and AUPRC.

import pickle
import numpy as np
import pandas as pd
import random
import math
from sklearn import metrics

def cal_AUROC(sorted_gene_list, positive_set):
	total = len(sorted_gene_list)
	actual_positive = positive_set & set(sorted_gene_list)
	# x: FPR = FP / Actual Negative
	# y: TPR = TP / Actual Positive
	AP = len(actual_positive)
	AN = total - AP
	FP = 0
	TP = 0
	coordinates = []
	for gene in sorted_gene_list:
		if gene in actual_positive:
			TP += 1
		else:
			FP += 1
		coordinates.append((FP/AN, TP/AP))
	FPR = [x[0] for x in coordinates]
	TPR = [x[1] for x in coordinates]	
	AUC = metrics.auc(FPR, TPR)
	return AUC

def cal_AUPRC(sorted_gene_list, positive_set):
	total = len(sorted_gene_list)
	actual_positive = positive_set & set(sorted_gene_list)
	# x: recall = (TP / Actual Positive)
	# y: precision = (TP / Predicted Positive)
	AP = len(actual_positive)
	PP = 0
	TP = 0
	coordinates = []
	for gene in sorted_gene_list:
		PP += 1
		if gene in actual_positive:
			TP += 1
		coordinates.append((TP/AP, TP/PP))
	recall = [x[0] for x in coordinates]
	precision = [x[1] for x in coordinates]
	AUC = metrics.auc(recall, precision)
	return AUC

def cal_AUPRC_percent(sorted_gene_list, positive_set, perc=1):
	k = math.ceil(len(sorted_gene_list) / 100 * perc)
	sorted_gene_list = sorted_gene_list[:k]
	total = len(sorted_gene_list)
	actual_positive = positive_set & set(sorted_gene_list)
	if len(actual_positive) == 0:
		return np.nan
	# x: recall = (TP / Actual Positive)
	# y: precision = (TP / Predicted Positive)
	AP = len(actual_positive)
	PP = 0
	TP = 0
	coordinates = []
	for gene in sorted_gene_list:
		PP += 1
		if gene in actual_positive:
			TP += 1
		coordinates.append((TP/AP, TP/PP))
	recall = [x[0] for x in coordinates]
	precision = [x[1] for x in coordinates]
	AUC = metrics.auc(recall, precision)
	return AUC

# GOBP_level: GOBP -> level
with open("GOBP_level.pkl", "rb") as INPUT:
	GOBP_level = pickle.load(INPUT)

# GOBP_BPancestors: GOBP -> Ancestor GOBP terms
with open("GOBP_BP_ancestors.pkl", "rb") as INPUT:
	GOBP_BPancestors = pickle.load(INPUT)

ensID_BPterms = {}
with open("EnsID_GOid_GOterm_GOdomain.BP_only.tsv") as INPUT:
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
with open("HumanRefProteome_UniProtACC_EnsGeneID.tsv") as INPUT:
	for line in INPUT:
		ensID = line.strip().split("\t")[1]
		if ensID != "N/A":
			ensID_list = ensID.split(",")
			reference_proteome.update(ensID_list)

df = pd.read_csv("CellxGene/ratio_of_n_and_n_cells_cell_type.gene.tissue-celltype.tsv", sep="\t", header=0, index_col=0)

# all proteins considered, so the intersection between dataframe and the human reference proteome
domain = list(set(df.index) & reference_proteome)
print("domain:", len(domain))

gene_gene2_corr = {}
na_list = []
for gene in domain:
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

print("gene pairs with corr value:", len(gene_gene2_corr))
print("gene pairs with NA value:", len(na_list))

gene_sorted_corr = dict(sorted(gene_gene2_corr.items(), key=lambda item: item[1], reverse=True))
del gene_gene2_corr

positive_set = set() # "gene2" that shares at least 1 GOBP (level>2) with "gene"
total = []

for (gene1, gene2) in list(gene_sorted_corr) + na_list:

	# If a gene doesn't have any GOBP (level > 2) annotations...
	# Should we include mark this gene as negative or exclude from the calculation?
	if gene1 not in ensID_BPterms or gene2 not in ensID_BPterms:
		continue
	else:
		total.append((gene1,gene2))

	gene1_BP_set = ensID_BPterms[gene1]
	gene2_BP_set = ensID_BPterms[gene2]
	count = len(gene1_BP_set & gene2_BP_set)
	if count > 0:
		positive_set.add((gene1, gene2))

print("positive gene pairs:", len(positive_set))
print("total gene pairs considered:", len(total))

if len(positive_set) > 0:
	AUROC = cal_AUROC(total, positive_set)
	AUPRC = cal_AUPRC(total, positive_set)
	AUPRC_1p = cal_AUPRC_percent(total, positive_set, 1)
else:
	AUROC = np.nan
	AUPRC = np.nan
	AUPRC_1p = np.nan

print(f"AUROC: {AUROC}, AUPRC: {AUPRC}, AUPRC_1%: {AUPRC_1p}")
