import pandas as pd
import numpy as np
from collections import defaultdict
import networkx as nx
from sklearn.model_selection import GroupShuffleSplit


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
            lambda row: frozenset([row['Gene_A_complex'], row['Gene_B_complex']]), axis=1
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
            lambda row: complex_to_supergroup.get(row['Gene_A_complex'], -1), axis=1
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

