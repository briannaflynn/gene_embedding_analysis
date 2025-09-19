import pandas as pd
import numpy as np
from collections import defaultdict
import networkx as nx
from sklearn.model_selection import GroupShuffleSplit

def assign_same_complex_multi(df):
    """
    Assigns Same_Complex_Multi = 1 for gene pairs that co-occur in at least one complex.
    Adds a unique gene pair ID and preserves all original rows.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns:
            - 'Gene_A', 'Gene_B'
            - 'complex_id_A', 'complex_id_B'

    Returns
    -------
    pd.DataFrame
        Original DataFrame with:
            - 'pair_id': unique unordered gene pair ID
            - 'pairwise_label': 1 if complex_id_A == complex_id_B else 0
            - 'Same_Complex_Multi': 1 if the pair co-occurred in any complex, else 0
    """
    df = df.copy()

    # Canonicalize gene pair order
    df['Gene_1'] = df[['Gene_A', 'Gene_B']].min(axis=1)
    df['Gene_2'] = df[['Gene_A', 'Gene_B']].max(axis=1)

    display(df)

    # Create a unique pair ID
    df['pair_id'] = df['Gene_1'] + '_' + df['Gene_2']

    #display(df)
    display(df[['Gene_1', 'Gene_2', 'Gene_A', 'Gene_B', 'pair_id']])

    # Label whether this row is a within-complex match
    df['pairwise_label'] = (df['complex_id_A'] == df['complex_id_B']).astype(int)
    display(df)

    # Flag any pair ID with at least one match
    pair_to_label = df.groupby('pair_id')['pairwise_label'].sum() >= 1
    df['Same_Complex_Multi'] = df['pair_id'].map(pair_to_label).astype(int)
    display(df[['Gene_1', 'Gene_2', 'Gene_A', 'Gene_B', 'pair_id', 'pairwise_label', 'Same_Complex_Multi']])

    return df


class ProteinComplexGroupSplitter:
    def __init__(self, positive_label=1, negative_ratio=3, random_state=42):
        self.positive_label = positive_label
        self.negative_ratio = negative_ratio
        self.random_state = random_state
        self.supergroup_map = {}
        self.positive_df = None
        self.negative_df = None
        self.full_df = None

    def _extract_positive_pairs(self, df):
        pos_df = df[df['Same_Complex'] == self.positive_label].copy()
        pos_df['complex_pair'] = pos_df.apply(
            lambda row: frozenset([row['complex_id_A'], row['complex_id_B']]), axis=1
        )
        return pos_df

    def _build_supergroups(self, pos_df):
        g = nx.Graph()
        for cpair in pos_df['complex_pair']:
            g.add_nodes_from(cpair)
            for a in cpair:
                for b in cpair:
                    if a != b:
                        g.add_edge(a, b)

        supergroups = list(nx.connected_components(g))
        complex_to_supergroup = {}
        for i, sg in enumerate(supergroups):
            for complex_id in sg:
                complex_to_supergroup[complex_id] = i

        pos_df['group_id'] = pos_df.apply(
            lambda row: complex_to_supergroup.get(row['complex_id_A'], -1), axis=1
        )
        return pos_df, complex_to_supergroup

    def _assign_negative_groups(self, df, pos_group_ids):
        neg_df = df[df['Same_Complex'] != self.positive_label].copy()
        sample_size = min(len(neg_df), self.negative_ratio * len(self.positive_df))
        neg_df = neg_df.sample(n=sample_size, random_state=self.random_state).copy()
        neg_df['group_id'] = np.random.choice(pos_group_ids, size=len(neg_df), replace=True)
        return neg_df

    def prepare(self, df):
        self.positive_df = self._extract_positive_pairs(df)
        self.positive_df, self.supergroup_map = self._build_supergroups(self.positive_df)
        pos_group_ids = self.positive_df['group_id'].unique()
        self.negative_df = self._assign_negative_groups(df, pos_group_ids)
        self.full_df = pd.concat([self.positive_df, self.negative_df], ignore_index=True)
        return self.full_df

    def group_split(self, test_size=0.25):
        if self.full_df is None:
            raise ValueError("Run `prepare(df)` first to generate group assignments.")

        X = self.full_df.drop(columns=['Same_Complex'])
        y = self.full_df['Same_Complex']
        groups = self.full_df['group_id']

        gss = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=self.random_state)
        for train_idx, test_idx in gss.split(X, y, groups):
            return (
                X.iloc[train_idx], X.iloc[test_idx],
                y.iloc[train_idx], y.iloc[test_idx],
                groups.iloc[train_idx], groups.iloc[test_idx]
            )


def stratified_group_split_balanced_per_fold(df, label_col='Same_Complex_Multi', group_col='pair_id',
                                             n_splits=5, random_state=42):
    """
    Stratified group-based CV with balanced class sampling per fold (equal pos/neg in train/test).

    Parameters
    ----------
    df : pd.DataFrame
        Dataset with row-level data.
    label_col : str
        Column name for binary label.
    group_col : str
        Group identifier column (e.g. 'pair_id').
    n_splits : int
        Number of folds.
    random_state : int
        Random seed.

    Returns
    -------
    List of (train_idx, test_idx) tuples with class-balanced samples per fold.
    """
    from sklearn.model_selection import StratifiedKFold
    import numpy as np

    rng = np.random.RandomState(random_state)

    # Create group-level label
    group_df = df.groupby(group_col)[label_col].max().reset_index()

    # Stratified KFold by group-level labels
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    splits = []

    for fold_idx, (train_group_idx, test_group_idx) in enumerate(skf.split(group_df, group_df[label_col])):
        train_groups = group_df.iloc[train_group_idx]
        test_groups = group_df.iloc[test_group_idx]

        # Expand to full row-level DataFrames
        train_rows = df[df[group_col].isin(train_groups[group_col])]
        test_rows = df[df[group_col].isin(test_groups[group_col])]

        # Split by label
        train_pos = train_rows[train_rows[label_col] == 1]
        train_neg = train_rows[train_rows[label_col] == 0]
        test_pos = test_rows[test_rows[label_col] == 1]
        test_neg = test_rows[test_rows[label_col] == 0]

        # Downsample majority class
        n_train = min(len(train_pos), len(train_neg))
        n_test = min(len(test_pos), len(test_neg))

        train_pos_sampled = train_pos.sample(n=n_train, random_state=random_state + fold_idx)
        train_neg_sampled = train_neg.sample(n=n_train, random_state=random_state + fold_idx + 1)
        test_pos_sampled = test_pos.sample(n=n_test, random_state=random_state + fold_idx + 2)
        test_neg_sampled = test_neg.sample(n=n_test, random_state=random_state + fold_idx + 3)

        # Combine and index
        train_idx = pd.concat([train_pos_sampled, train_neg_sampled]).index.tolist()
        test_idx = pd.concat([test_pos_sampled, test_neg_sampled]).index.tolist()

        splits.append((train_idx, test_idx))

    return splits


folds = stratified_group_split_balanced_per_fold(
    df=prepared_df,
    label_col='Same_Complex',
    group_col='group_id',
    n_splits=5,
    random_state=42
)

for i, (train_idx, test_idx) in enumerate(folds):

    tr = prepared_df.loc[train_idx]
    tr['Test'] = False
    ts = prepared_df.loc[test_idx]
    ts['Test'] = True

    # data = pd.concat([tr, ts])
    # data.drop(columns = ['complex_pair'], inplace=True)
    # data.to_csv(f'balanced_group_stratified_fold_{i}.csv',index=False)

    y_train = prepared_df.loc[train_idx, 'Same_Complex']
    y_test = prepared_df.loc[test_idx, 'Same_Complex']
    print(f"Fold {i+1}:")
    print(f"  Train: {len(train_idx)} rows ({y_train.sum()} pos, {len(train_idx) - y_train.sum()} neg)")
    print(f"  Test:  {len(test_idx)} rows ({y_test.sum()} pos, {len(test_idx) - y_test.sum()} neg)")


