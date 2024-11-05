#!/usr/bin/python

import pandas as pd
import numpy as np
from scipy.spatial.distance import cosine

def ensembl_to_gene_name(ensembl_id):
    """
    Converts a single Ensembl ID to the corresponding gene name.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        'query': f'{ensembl_id} AND organism_id:9606',  # Filtering for human (organism ID 9606)
        'fields': 'gene_names',  # Fetch gene name
        'format': 'tsv',         # Requesting tab-separated values format
        'size': 1                # Limit results to one entry
    }
    response = requests.get(url, params=params)
    if response.status_code == 200:
        lines = response.text.splitlines()
        #print(lines)
        #print(len(lines))
         # Check if we have data in the second line
        if len(lines) > 1:
            gene_name_data = lines[1].split('\t')
            if len(gene_name_data) > 0:  # Check if there’s at least one element in gene_name_data
                gene_name = gene_name_data[0]
                return gene_name
    return None

def convert_ensembl_ids_and_compute_similarity(df, embeddings_dict):
    """
    Converts pairs of Ensembl IDs to gene names, computes cosine similarity
    between corresponding embedding vectors, and returns a merged DataFrame.
    
    Args:
        df (pd.DataFrame): DataFrame with 'ensembl_id_pair' and 'correlation' columns.
        embeddings_dict (dict): Dictionary where keys are gene names, values are embedding vectors.
        
    Returns:
        pd.DataFrame: A DataFrame with Ensembl ID pairs, gene name pairs, Pearson correlations,
                      cosine similarities, and any additional columns from the original DataFrame.
    """
    # Initialize lists to store results
    ensembl_id_1, ensembl_id_2 = [], []
    gene_name_1, gene_name_2 = [], []
    cosine_similarities = []

    # Process each row in the DataFrame
    for idx, row in df.iterrows():
        ensembl_pair = row['ensembl_id_pair']
        ensembl_id1, ensembl_id2 = ensembl_pair.strip("()").replace("'", "").split(", ")

        print(ensembl_id1, ensembl_id2)
        
        # Convert Ensembl IDs to gene names
        gene1 = ensembl_to_gene_name(ensembl_id1)
        gene2 = ensembl_to_gene_name(ensembl_id2)

        # Check if gene names are missing
        if gene1 is None or gene2 is None:
            print(f"No associated gene name found for Ensembl ID pair: {ensembl_id1}, {ensembl_id2}")

        print(gene1, gene2)
        
        # Append Ensembl IDs and gene names
        ensembl_id_1.append(ensembl_id1)
        ensembl_id_2.append(ensembl_id2)
        gene_name_1.append(gene1)
        gene_name_2.append(gene2)
        
        # Compute cosine similarity if both gene names are found in embeddings_dict
        if gene1 in embeddings_dict and gene2 in embeddings_dict:
            embedding1 = embeddings_dict[gene1]
            embedding2 = embeddings_dict[gene2]
            cos_sim = 1 - cosine(embedding1, embedding2)  # Cosine similarity is 1 - cosine distance
        else:
            cos_sim = np.nan  # Use NaN if one or both gene names are missing in embeddings_dict
        cosine_similarities.append(cos_sim)

    # Add columns to the DataFrame
    result_df = df.copy()
    result_df['ensembl_id_1'] = ensembl_id_1
    result_df['ensembl_id_2'] = ensembl_id_2
    result_df['gene_name_1'] = gene_name_1
    result_df['gene_name_2'] = gene_name_2
    result_df['cosine_similarity'] = cosine_similarities

    return result_df
