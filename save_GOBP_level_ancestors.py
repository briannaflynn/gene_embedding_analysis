#!/usr/bin/env python3
#Date: September 4, 2024; September 30, 2024
#Author: Muyoung Lee
#Description: Print [GOBP term ID, Level, Ancestor BP terms].
#Usage: (scanpy) [THIS SCRIPT]

#Ref: https://github.com/tanghaibao/goatools/blob/main/notebooks/parents_and_ancestors.ipynb
#Ref2: https://github.com/tanghaibao/goatools/blob/main/notebooks/report_depth_level.ipynb

import pickle

from goatools.obo_parser import GODag
godag = GODag("go-basic.obo.1", optional_attrs={"relationsihp"})

from goatools.gosubdag.gosubdag import GoSubDag
#optional_relationships = {'regulates', 'negatively_regulates', 'positively_regulates'}

BP_level = {}
BP_ancestors = {}

for GOID in godag:
	if godag[GOID].namespace != "biological_process":
		continue
	BP_level[GOID] = godag[GOID].level
	if godag[GOID].name == "biological_process":
		continue
	ancestors = list(GoSubDag([GOID], godag, prt=None).rcntobj.go2ancestors[GOID])
	ancestors_BP = []
	for go in ancestors:
		if godag[go].namespace == "biological_process":
			ancestors_BP.append(go)
	BP_ancestors[GOID] = ancestors_BP

with open("GOBP_level.pkl", "wb") as OUTPUT1:
	pickle.dump(BP_level, OUTPUT1)

with open("GOBP_BP_ancestors.pkl", "wb") as OUTPUT2:
	pickle.dump(BP_ancestors, OUTPUT2)
