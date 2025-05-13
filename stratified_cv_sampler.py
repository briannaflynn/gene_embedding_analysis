import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold

def stratified_downsample_cv(df, label_col='Label', n_splits=5, random_state=42, return_idx = False):
    """
    Performs stratified k-fold cross-validation with downsampling of the majority class.

    BALANCED TRAINING SET, BALANCED TEST SET

    Parameters:
        df (pd.DataFrame): The input dataset containing features and a binary label column.
        label_col (str): The name of the label column.
        n_splits (int): Number of folds.
        random_state (int): Random seed for reproducibility.

    Returns:
        List of (train_df, test_df) tuples, each representing a fold.
    """
    df_positive = df[df[label_col] == 1].copy()
    df_negative = df[df[label_col] == 0].copy()

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    folds = list(skf.split(df_positive, df_positive[label_col]))

    fold_data = []

    if return_idx:

        for fold_idx, (train_pos_idx, test_pos_idx) in enumerate(folds):
            train_pos = df_positive.iloc[train_pos_idx]
            test_pos = df_positive.iloc[test_pos_idx]
    
            # Downsample negatives
            train_neg = df_negative.sample(n=len(train_pos), random_state=random_state + fold_idx)
            test_neg = df_negative.drop(train_neg.index, errors='ignore').sample(
                n=len(test_pos), random_state=random_state + 10 + fold_idx
            )
    
            # Get indices
            train_idx = pd.concat([train_pos, train_neg]).index.to_list()
            test_idx = pd.concat([test_pos, test_neg]).index.to_list()
    
            fold_data.append((train_idx, test_idx))

    else:
        
        for fold_idx, (train_pos_idx, test_pos_idx) in enumerate(folds):
            train_pos = df_positive.iloc[train_pos_idx]
            test_pos = df_positive.iloc[test_pos_idx]
    
            # Downsample negatives
            train_neg = df_negative.sample(n=len(train_pos), random_state=random_state + fold_idx)
            test_neg = df_negative.drop(train_neg.index, errors='ignore').sample(
                n=len(test_pos), random_state=random_state + 10 + fold_idx
            )
    
            train_df = pd.concat([train_pos, train_neg]).sample(frac=1.0, random_state=random_state + 20 + fold_idx)
            test_df = pd.concat([test_pos, test_neg]).sample(frac=1.0, random_state=random_state + 30 + fold_idx)
    
            fold_data.append((train_df.reset_index(drop=True), test_df.reset_index(drop=True)))

    return fold_data

def stratified_cv_full_neg_train(df, label_col='Label', n_splits=5, random_state=42, return_idx=False):
    """
    Performs stratified k-fold cross-validation:

    UNBALANCED TRAINING SET, BALANCED TEST SET
    - Training set includes all negative samples and positive samples in the training split.
    - Test set is balanced: equal number of positives and downsampled negatives.

    Parameters:
        df (pd.DataFrame): The input dataset.
        label_col (str): The name of the label column.
        n_splits (int): Number of folds.
        random_state (int): Random seed for reproducibility.

    Returns:
        List of (train_df, test_df) tuples.
    """
    df_positive = df[df[label_col] == 1].copy()
    df_negative = df[df[label_col] == 0].copy()

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    folds = list(skf.split(df_positive, df_positive[label_col]))

    fold_data = []

    if return_idx:
    
        for fold_idx, (train_pos_idx, test_pos_idx) in enumerate(folds):
            train_pos = df_positive.iloc[train_pos_idx]
            test_pos = df_positive.iloc[test_pos_idx]
    
            # Training set: all negatives + train positives
            train_idx = pd.concat([train_pos, df_negative]).index.to_list()
    
            # Test set: balanced test positives + sampled negatives
            test_neg = df_negative.sample(n=len(test_pos), random_state=random_state + 10 + fold_idx)
            test_idx = pd.concat([test_pos, test_neg]).index.to_list()
    
            fold_data.append((train_idx, test_idx))

    else:
    
        for fold_idx, (train_pos_idx, test_pos_idx) in enumerate(folds):
            # Train: all negatives + positive train split
            train_pos = df_positive.iloc[train_pos_idx]
            train_df = pd.concat([train_pos, df_negative]).sample(frac=1.0, random_state=random_state + fold_idx)
    
            # Test: positive test split + downsampled negatives
            test_pos = df_positive.iloc[test_pos_idx]
            test_neg = df_negative.sample(n=len(test_pos), random_state=random_state + 10 + fold_idx)
            test_df = pd.concat([test_pos, test_neg]).sample(frac=1.0, random_state=random_state + 20 + fold_idx)
    
            fold_data.append((train_df.reset_index(drop=True), test_df.reset_index(drop=True)))

    return fold_data

def generate_test_data(n_pos=10000, n_neg=2000000, n_features=5):
    np.random.seed(0)
    # Create features for positives and negatives
    pos_data = pd.DataFrame(np.random.randn(n_pos, n_features), columns=[f"feature_{i}" for i in range(n_features)])
    pos_data['Label'] = 1

    neg_data = pd.DataFrame(np.random.randn(n_neg, n_features), columns=[f"feature_{i}" for i in range(n_features)])
    neg_data['Label'] = 0

    full_df = pd.concat([pos_data, neg_data]).sample(frac=1.0).reset_index(drop=True)
    return full_df
