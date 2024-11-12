#!/usr/bin/env python3
#Date: Sept 3, Oct 4, Oct 14, Oct 23, Nov 11, 2024
#Author: Muyoung Lee
#Description:
#1. From the reference human proteome genes, gather genes belonging to the input feature table.
#   Exclude genes less than 30 non-zero values (observations).
#   Exclude genes which have only too broad GOBP terms.
#2. Calculate Spearman correlataions between all gene pairs.
#3. Sort gene pairs based on calculated correlations.
#4. If the two genes in one gene pair share at least one GOBP term whose level is larger than 2 (deeper than 2), mark this gene pair as positive. Consider ancestors GOBP terms as well.
#5. Based on steps 3 and 4, calculate AUROC and AUPRC.
#Usage: [THIS SCRIPT] [GOBP_term_size_threshold]

import sys
import pickle
import numpy as np
import pandas as pd
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
	actual_positive = positive_set & set(sorted_gene_list)
	if len(actual_positive) == 0:
		return np.nan
	else:
		return cal_AUPRC(sorted_gene_list, positive_set)

# GOBP_level: GOBP -> level
with open("GOBP/GOBP_level.pkl", "rb") as INPUT:
	GOBP_level = pickle.load(INPUT)

# GOBP_BPancestors: GOBP -> Ancestor GOBP terms
with open("GOBP/GOBP_BP_ancestors.pkl", "rb") as INPUT:
	GOBP_BPancestors = pickle.load(INPUT)

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

large_BPterms = set()
with open("GOBP/GOBP_size_reference_proteome_genes.tsv") as INPUT:
	for line in INPUT:
		GOBP, size = line.strip().split("\t")[:2]
		if int(size) >= int(sys.argv[1]):
			large_BPterms.add(GOBP)

print(f"GOBP term size threshold: {sys.argv[1]}", file=sys.stderr)

adj_ensID_BPterms = {}
for ensID in ensID_BPterms:
	new_set = ensID_BPterms[ensID] - large_BPterms
	if len(new_set) > 0:
		adj_ensID_BPterms[ensID] = new_set
del ensID_BPterms

df = pd.read_csv("CZ_CellxGene/ratio_of_n_and_n_cells_cell_type.gene.tissue-celltype.ver2.tsv", sep="\t", header=0, index_col=0)

# all proteins considered, so the intersection between dataframe and the human reference proteome
domain = list(set(df.index) & reference_proteome)

too_many_zero_genes = set()
for gene in domain:
	non_zero_count = (df.loc[gene] != 0).sum()
	if non_zero_count < 30:
		too_many_zero_genes.add(gene)

print("genes with less than 30 observations:", len(too_many_zero_genes), file=sys.stderr)
domain = set(domain) - too_many_zero_genes
print(f"gene with at least 30 observations: {len(domain)}", file=sys.stderr)

domain = domain & set(adj_ensID_BPterms)
print(f"gene with at least 30 observations and GOBP terms: {len(domain)}", file=sys.stderr)

gene_gene2_corr = {}
for gene in domain:
	for gene2 in domain:
		if gene2 == gene:
			continue
		gene_, gene2_ = sorted([gene, gene2])
		if (gene_, gene2_) in gene_gene2_corr:
			continue
		cor_value = df.loc[gene_].corr(df.loc[gene2_], method="spearman")
		gene_gene2_corr[(gene_, gene2_)] = cor_value

gene_sorted_corr = dict(sorted(gene_gene2_corr.items(), key=lambda item: item[1], reverse=True))
del gene_gene2_corr

positive_set = set() # Gene pairs that share at least 1 GOBP (level>2) term.
total = list(gene_sorted_corr)

for (gene1, gene2) in gene_sorted_corr:
	gene1_BP_set = adj_ensID_BPterms[gene1]
	gene2_BP_set = adj_ensID_BPterms[gene2]
	count = len(gene1_BP_set & gene2_BP_set)
	if count > 0:
		positive_set.add((gene1, gene2))

print(f"number of genes pairs considered: {len(total)}", file=sys.stderr)
print("positive gene pairs:", len(positive_set), file=sys.stderr)

for gene_pair in total:
	is_positive = 0
	if gene_pair in positive_set:
		is_positive = 1
	print(gene_pair, gene_sorted_corr[gene_pair], is_positive, sep="\t")

if len(positive_set) > 0:
	AUROC = cal_AUROC(total, positive_set)
	AUPRC = cal_AUPRC(total, positive_set)
	AUPRC_1p = cal_AUPRC_percent(total, positive_set, 1)
else:
	AUROC = np.nan
	AUPRC = np.nan
	AUPRC_1p = np.nan

print(f"AUROC: {AUROC:.3f}, AUPRC: {AUPRC:.3f}, AUPRC_1%: {AUPRC_1p:.3f}", file=sys.stderr)
