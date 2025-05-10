import pandas as pd

def create_embedding_feature_matrix(pair_df, embedding_dict):
    """
    Constructs a feature matrix for gene pairs using their embeddings.

    Parameters:
        pair_df (pd.DataFrame): Must contain 'Gene_A', 'Gene_B', 'Label'.
        embedding_dict (dict): gene_name -> embedding (list or np.array)

    Returns:
        pd.DataFrame: Rows indexed by 'GeneA_GeneB', columns are embedding dims + label
    """
    records = []
    for _, row in pair_df.iterrows():
        gene_a = row['Gene_A']
        gene_b = row['Gene_B']
        label = row['Label']

        if gene_a not in embedding_dict or gene_b not in embedding_dict:
            continue  # Skip if embeddings are missing

        emb_a = embedding_dict[gene_a]
        emb_b = embedding_dict[gene_b]

        record = {
            **{f'gene_a_{i}': val for i, val in enumerate(emb_a)},
            **{f'gene_b_{i}': val for i, val in enumerate(emb_b)},
            'Label': label
        }

        records.append((f"{gene_a}_{gene_b}", record))

    feature_df = pd.DataFrame.from_dict(dict(records), orient='index')
    feature_df.index.name = "Gene_Pair"

    return feature_df
