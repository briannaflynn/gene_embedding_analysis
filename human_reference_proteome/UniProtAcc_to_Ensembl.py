#!/usr/bin/env python3
#Date: Oct 1, 2024
#Author: Muyoung Lee
#Description: Print a table of UniProt-accession and EnsemblID for the reference human proteome
#Usage: [THIS SCRIPT]

acc_hgnc = {}
with open("uniprotkb_proteome_UP000005640_AND_revi_2024_09_30.tsv") as INPUT:
	INPUT.readline()
	for line in INPUT:
		words = line.strip("\n").split("\t")
		acc = words[0]
		hgnc = words[3]
		if acc != "" and hgnc != "":
			hgnc_list = []
			for i in hgnc.split(";"):
				if i != "":
					hgnc_list.append(i)
		acc_hgnc[acc] = hgnc_list

hgnc_ens = {}
with open("HGNC_Ensembl.txt") as INPUT:
	INPUT.readline()
	for line in INPUT:
		words = line.strip("\n").split("\t")
		hgnc, ensID = words[0], words[1]
		if ensID != "":
			hgnc_ens[hgnc] = ensID

for acc in acc_hgnc:
	ens_set = set()
	for hgnc in acc_hgnc[acc]:
		if hgnc in hgnc_ens:
			ens_set.add(hgnc_ens[hgnc])
	if len(ens_set) > 0:
		print(acc, ",".join(sorted(list(ens_set))), sep="\t")
	else:
		print(acc, "N/A", sep="\t")
