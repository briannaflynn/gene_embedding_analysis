#!/usr/bin/env python3
#Date: Dec 04, 2024
#Author: Brianna Flynn
#Description:
#This script is for processing different thresholds from 150-1500 to test how this effects "domain size", print this to a table and save each dictionary
import sys
import pickle
import numpy as np
import pandas as pd
import math

# python prep_domain.py 100 1500 2
gobp_term_size_threshold = int(sys.argv[1])
gobp_count_celltype_threshold = int(sys.argv[2])
gobp_level = int(sys.argv[3])

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
		elif GOBP_level[goid] <= gobp_level:
			continue
	
		if ensID not in ensID_BPterms:
			ensID_BPterms[ensID] = set([goid])
		else:
			ensID_BPterms[ensID].add(goid)

		if goid in GOBP_BPancestors:
			ancestors = []
			for GOBP in GOBP_BPancestors[goid]:
				if GOBP_level[GOBP] > gobp_level:
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

large_BPterms = set()
with open("BPterm_size_genes.tsv") as INPUT:
	for line in INPUT:
		GOBP, size = line.strip().split("\t")[:2]
		if int(size) >= gobp_term_size_threshold:
			large_BPterms.add(GOBP)

print(f"GOBP term size threshold: {gobp_term_size_threshold}", file=sys.stderr)

# Efficiently filter ensID_BPterms using dictionary comprehension
adj_ensID_BPterms = {
    ensID: terms - large_BPterms
    for ensID, terms in ensID_BPterms.items()
    if len(terms - large_BPterms) > 0
}

df = pd.read_csv("ratio_of_n_and_n_cells_cell_type.gene.tissue-celltype.ver2.tsv", sep="\t", header=0, index_col=0)

# all proteins considered, so the intersection between dataframe and the human reference proteome
domain = list(set(df.index) & reference_proteome)

domain_df = pd.DataFrame()
domain_df['Domain']=domain
domain_df.to_csv('all_proteins_considered.intersection_df_human-ref-proteome.csv', index=False)

# # applying the count filtering now
valid_genes = df.loc[domain].apply(lambda x: (x != 0).sum(), axis=1) >= gobp_count_celltype_threshold
domain = set(df.loc[domain][valid_genes].index)

# update domain after filtering from gobp term size threshold to only include genes intersecting adj_ensID_BPterms
domain = domain & set(adj_ensID_BPterms)

print(f'After filtering based on a gobp term size threshold of {gobp_term_size_threshold} with at least {gobp_count_celltype_threshold} observations: {len(domain)}')

export = pd.DataFrame()
export['Domain'] = list(domain)
export.to_csv(f'all_proteins_considered.intersection_df_human-ref-proteome_size-{gobp_term_size_threshold}_level{gobp_level}_obs-{gobp_count_celltype_threshold}.csv', index=False)
print(export.shape)