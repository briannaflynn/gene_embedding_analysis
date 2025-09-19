import sys
import os
import pandas as pd
import numpy as np
import pytest

# Add parent directory (where protein_complex.py lives) to sys.path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from protein_complex import (
    assign_same_complex_multi,
    ProteinComplexGroupSplitter,
    stratified_group_split_balanced_per_fold,
)


def test_assign_same_complex_multi_basic():
    df = pd.DataFrame({
        "Gene_A": ["A", "A", "B"],
        "Gene_B": ["B", "C", "C"],
        "complex_id_A": [1, 1, 2],
        "complex_id_B": [1, 2, 2],
    })

    result = assign_same_complex_multi(df)

    # Check added columns
    assert "pair_id" in result.columns
    assert "pairwise_label" in result.columns
    assert "Same_Complex_Multi" in result.columns

    # Check canonicalization
    assert all(result["pair_id"].str.contains("_"))

    # "A_B" row where both complex IDs = 1 should be labeled positive
    ab_row = result[result["pair_id"] == "A_B"]
    assert ab_row["pairwise_label"].iloc[0] == 1
    assert ab_row["Same_Complex_Multi"].iloc[0] == 1


def make_dummy_df():
    """Utility to build a small df with Same_Complex labels."""
    return pd.DataFrame({
        "Gene_A": ["A", "B", "C", "D"],
        "Gene_B": ["B", "C", "D", "E"],
        "complex_id_A": [1, 1, 2, 3],
        "complex_id_B": [1, 2, 2, 4],
        "Same_Complex": [1, 0, 1, 0],
    })


def test_group_splitter_prepare_and_split():
    df = make_dummy_df()
    splitter = ProteinComplexGroupSplitter(positive_label=1, negative_ratio=1, random_state=0)

    prepared = splitter.prepare(df)

    # Must contain positive and negative examples
    assert prepared["Same_Complex"].nunique() == 2
    assert "group_id" in prepared.columns

    # Group split should work
    X_train, X_test, y_train, y_test, g_train, g_test = splitter.group_split(test_size=0.5)
    assert not X_train.empty
    assert not X_test.empty
    # Train/test groups should not overlap
    assert set(g_train).isdisjoint(set(g_test))


def test_group_splitter_errors_without_prepare():
    splitter = ProteinComplexGroupSplitter()
    with pytest.raises(ValueError):
        splitter.group_split()


def test_stratified_group_split_balanced_per_fold():
    df = pd.DataFrame({
        "pair_id": ["A_B", "A_B", "C_D", "C_D", "E_F", "E_F", "G_H", "G_H"],
        "Same_Complex_Multi": [1, 1, 1, 1, 0, 0, 0, 0],
    })

    splits = stratified_group_split_balanced_per_fold(df, label_col="Same_Complex_Multi", group_col="pair_id", n_splits=2)

    assert len(splits) == 2  # 2 folds
    for train_idx, test_idx in splits:
        # Check indices are valid
        assert all(i in df.index for i in train_idx)
        assert all(i in df.index for i in test_idx)

        # Class balance should be equal
        train_labels = df.loc[train_idx, "Same_Complex_Multi"]
        test_labels = df.loc[test_idx, "Same_Complex_Multi"]
        assert train_labels.sum() == len(train_labels) // 2
        assert test_labels.sum() == len(test_labels) // 2

