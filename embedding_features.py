import numpy as np
import pandas as pd

def create_embedding_feature_matrix(pair_df, embedding_dict, gene_columns=('Gene_A', 'Gene_B'), label_column='Same_Complex'):
    """
    Efficiently constructs a feature matrix for gene pairs using their embeddings.
    """
    gene_a_col, gene_b_col = gene_columns
    embedding_dim = len(next(iter(embedding_dict.values())))
    
    valid_pairs = []
    features = []
    labels = []

    for row in pair_df.itertuples(index=False):
        gene_a = getattr(row, gene_a_col)
        gene_b = getattr(row, gene_b_col)
        label = getattr(row, label_column)

        emb_a = embedding_dict.get(gene_a)
        emb_b = embedding_dict.get(gene_b)

        if emb_a is None or emb_b is None:
            continue

        combined = np.concatenate([emb_a, emb_b])
        features.append(combined)
        labels.append(label)
        valid_pairs.append(f"{gene_a}_{gene_b}")

    feature_array = np.array(features)
    feature_columns = [f"gene_a_{i}" for i in range(embedding_dim)] + [f"gene_b_{i}" for i in range(embedding_dim)]
    
    feature_df = pd.DataFrame(feature_array, columns=feature_columns, index=valid_pairs)
    feature_df["Label"] = labels
    feature_df.index.name = "Gene_Pair"

    return feature_df
