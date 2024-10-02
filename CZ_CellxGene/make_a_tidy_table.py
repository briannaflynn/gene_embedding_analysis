#!/usr/bin/env python3
#Date: September 10, 2024
#Author: Muyoung Lee
#Description:
# 0. Find human and healthy samples.
# 1. Calculate the number of cells for specific cell type & tissue
# 2. Calculate the number of cells with non-zero expression for each gene in a cell type & tissue pair.
#Usage: [THIS SCRIPT]

import csv

gene_tissue_cell_whole = {}
gene_tissue_cell_exp = {}
all_tissue_celltypes = {}

'''
0: index
1: tissue_ontology_term_id
2: organism_ontology_term_id
3: tissue_original_ontology_term_id
4: dataset_id
5: disease_ontology_term_id
6: self_reported_ethnicity_ontology_term_id
7: sex_ontology_term_id
8: publication_citation
9: gene_ontology_term_id
10: cell_type_ontology_term_id
11: n_cells_cell_type
12: n_cells_tissue
13: n
14: me
15: pc
16: tpc
17: tissue_name
18: organism_name
19: gene_name
20: cell_type_name
21: disease_name
22: self_reported_ethnicity_name
23: sex_name
24: tissue_original_name
'''

# https://cellxgene.cziscience.com/docs/04__Analyze%20Public%20Data/4_2__Gene%20Expression%20Documentation/4_2_6__Gene%20Expression%20Source%20Data
with open("expression-summary-full-03-11-24.csv") as INPUT:
	INPUT.readline()
	for line in csv.reader(INPUT):
		if line[2] != "NCBITaxon:9606" or line[21] != "normal":
			continue
		gene = line[9]
		n_cells_cell_type = int(line[11])
		n = int(line[13])
		tissue = line[17]
		cell_type = line[20]
		if gene not in gene_tissue_cell_whole:
			gene_tissue_cell_whole[gene] = {}
			gene_tissue_cell_exp[gene] = {}
		if tissue not in gene_tissue_cell_whole[gene]:
			gene_tissue_cell_whole[gene][tissue] = {}
			gene_tissue_cell_exp[gene][tissue] = {}
		if cell_type not in gene_tissue_cell_whole[gene][tissue]:
			gene_tissue_cell_whole[gene][tissue][cell_type] = 0
			gene_tissue_cell_exp[gene][tissue][cell_type] = 0
		gene_tissue_cell_whole[gene][tissue][cell_type] += n_cells_cell_type
		gene_tissue_cell_exp[gene][tissue][cell_type] += n

		if tissue not in all_tissue_celltypes:
			all_tissue_celltypes[tissue] = set()
		all_tissue_celltypes[tissue].add(cell_type)

header = "gene\t"
for tissue in sorted(all_tissue_celltypes):
	for celltype in sorted(all_tissue_celltypes[tissue]):
		pwd = f"{tissue.replace(' ', '_')}@{celltype.replace(' ', '_')}"
		header += f"{pwd}\t"
print(header.rstrip())

for gene in sorted(gene_tissue_cell_whole):
	outstr = f"{gene}\t"
	for tissue in sorted(all_tissue_celltypes):
		for celltype in sorted(all_tissue_celltypes[tissue]):
			if tissue in gene_tissue_cell_whole[gene] and celltype in gene_tissue_cell_whole[gene][tissue]:
				if gene_tissue_cell_whole[gene][tissue][celltype] != 0:
					value = gene_tissue_cell_exp[gene][tissue][celltype] / gene_tissue_cell_whole[gene][tissue][celltype]
				else:
					value = "N/A"
			else:
				value = "N/A"
			outstr += f"{value}\t"
	print(outstr.rstrip())
