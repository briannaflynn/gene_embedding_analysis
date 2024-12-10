#!/usr/bin/env python3
#Date: Dec 09, 2024
#Author: Brianna Flynn
#Description:
#Processes gene pairs from a dataframe and identifies pairs that share Gene Ontology Biological Process (GOBP) terms.
import pandas as pd
import pickle
from itertools import combinations
import sys 

def process_gene_pairs_from_dataframe(
    gobp_level_threshold, 
    gobp_size_threshold, 
    df, 
    ensID_GO_file, 
    bpterm_size_file, 
    gobp_level_file="GOBP_level.pkl", 
    gobp_ancestors_file="GOBP_BP_ancestors.pkl"
):
    """
    
    Parameters:
        gobp_level_threshold (int): Minimum GOBP level to consider.
        gobp_size_threshold (int): Maximum size of GOBP terms to exclude.
        df (pd.DataFrame): DataFrame where rows are gene IDs (index) and columns are features (not used in this function).
        ensID_GO_file (str): Path to file mapping Ensembl IDs to GO terms.
        bpterm_size_file (str): Path to file containing GOBP term sizes.
        gobp_level_file (str): Path to pickle file with GOBP levels.
        gobp_ancestors_file (str): Path to pickle file with GOBP ancestors.

    Returns:
        pd.DataFrame: A DataFrame with columns:
            - `gene_1`: First gene in the pair.
            - `gene_2`: Second gene in the pair.
            - `GOBP_value`: Binary value indicating if the gene pair shares a GOBP term (1) or not (0).
    
    Notes:
        - The function assumes the input dataframe `df` has genes as the index and columns can be ignored.
        - The gene pairs are generated using combinations of all unique rows (genes) in the dataframe.
        - Sharing of GOBP terms is determined based on the filtered terms after applying size and level thresholds.
    """
    # Load GOBP levels and ancestors
    with open(gobp_level_file, "rb") as f:
        GOBP_level = pickle.load(f)
    with open(gobp_ancestors_file, "rb") as f:
        GOBP_BPancestors = pickle.load(f)

    # Pre-filter GOBP levels
    valid_gobp = {goid: level for goid, level in GOBP_level.items() if level > gobp_level_threshold}
    print('valid gobp')
    #print(valid_gobp)
    #sys.exit()

    # Create mapping of Ensembl IDs to relevant GOBP terms
    ensID_BPterms = {}
    with open(ensID_GO_file) as f:
        for line in f:
            ensID, goid = line.strip().split("\t")[:2]
            if goid in valid_gobp:
                ensID_BPterms.setdefault(ensID, set()).add(goid)
                if goid in GOBP_BPancestors:
                    ancestors = {anc for anc in GOBP_BPancestors[goid] if anc in valid_gobp}
                    ensID_BPterms[ensID].update(ancestors)
    print('ens id bp terms')
    #print(ensID_BPterms)
    #print(dict(list(ensID_BPterms.items())[-10:]))
    #sys.exit()
    # Filter GOBP terms based on size
    large_BPterms = set()
    with open(bpterm_size_file) as f:
        for line in f:
            GOBP, size = line.strip().split("\t")[:2]
            if int(size) >= gobp_size_threshold:
                large_BPterms.add(GOBP)
    print('large bp terms')
    #print(list(large_BPterms)[-10:])
    #sys.exit()

    adj_ensID_BPterms = {
        ensID: terms - large_BPterms
        for ensID, terms in ensID_BPterms.items()
        if terms - large_BPterms
    }
    print('adj ens id bp terms')

    print('Precompute GOBP terms for each gene')
    gene_to_bp_terms = {gene: adj_ensID_BPterms.get(gene, set()) for gene in df.index}
    #print('precomputed gobp terms for each gene', gene_to_bp_terms)
    print('Generate all possible gene pairs')
    all_pairs = list(combinations(df.index, 2))
    #print('all pairs', all_pairs)
    print('Calculate positive set')
    positive_set = {
        (gene1, gene2)
        for gene1, gene2 in all_pairs
        if gene_to_bp_terms[gene1] & gene_to_bp_terms[gene2]
    }
    #print(positive_set)
    # Create a DataFrame directly
    data = {
        "gene_1": [pair[0] for pair in all_pairs],
        "gene_2": [pair[1] for pair in all_pairs],
        "GOBP_value": [
            1 if pair in positive_set else 0
            for pair in all_pairs
        ]
    }
    return pd.DataFrame(data)


df = pd.read_csv("ratio_of_n_and_n_cells_cell_type.gene.tissue-celltype.ver2.tsv", sep="\t", header=0, index_col=0)
domain = pd.read_csv("all_proteins_considered.intersection_df_human-ref-proteome.csv")
domain = domain['Domain'].to_list()
debug_domain = domain[:10]
df = df.loc[domain]

result = process_gene_pairs_from_dataframe(
    gobp_level_threshold=2,
    gobp_size_threshold=10000,
    df=df,
    ensID_GO_file="EnsID_GOid_GOterm_GOdomain.BP_only.tsv",
    bpterm_size_file="BPterm_size_genes.tsv",
    gobp_level_file="GOBP_level.pkl",
    gobp_ancestors_file="GOBP_BP_ancestors.pkl"
)

print(result)
print(result['GOBP_value'].sum())

result.to_pickle("gobp_size10000_level2.pkl")

# Total rows: 202779591
# LEVEL 2
# gobp 100 = 1654358 positive
# gobp 1000 = 18843112 positive
