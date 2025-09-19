# Protein Complex Pair Processing and Splitting

This module provides utilities for preparing datasets of gene/protein pairs, labeling them by complex co-membership, and performing group-aware stratified cross-validation.

---

## `assign_same_complex_multi(df)`

**Purpose**  
Annotates gene pairs with whether they belong to the same protein complex, adds canonicalized identifiers, and creates a pair-level label that accounts for multiple observations.

### Input Format
A `pandas.DataFrame` with row-level comparisons of two genes. Must contain:

- `Gene_A` (`str`): Identifier for the first gene.  
- `Gene_B` (`str`): Identifier for the second gene.  
- `complex_id_A` (`str` or `int`): Complex identifier for `Gene_A`.  
- `complex_id_B` (`str` or `int`): Complex identifier for `Gene_B`.

### Output Format
Original DataFrame with the following additional columns:

- `Gene_1`, `Gene_2`: Canonicalized ordering of the gene pair (`min`, `max`).  
- `pair_id`: Unique string ID of the form `"Gene_1_Gene_2"`.  
- `pairwise_label`: Row-level indicator (`1` if `complex_id_A == complex_id_B`, else `0`).  
- `Same_Complex_Multi`: Pair-level indicator (`1` if any row for that pair has `pairwise_label = 1`, else `0`).

### Example
| Gene_A | Gene_B | complex_id_A | complex_id_B | Gene_1 | Gene_2 | pair_id     | pairwise_label | Same_Complex_Multi |
|--------|--------|--------------|--------------|--------|--------|-------------|----------------|---------------------|
| SGCB   | SGCD   | C1           | C1           | SGCB   | SGCD   | SGCB_SGCD   | 1              | 1                   |
| SGCB   | SGCD   | C1           | C2           | SGCB   | SGCD   | SGCB_SGCD   | 0              | 1                   |
| PSMB1  | PSMA4  | C3           | C4           | PSMA4  | PSMB1  | PSMA4_PSMB1 | 0              | 0                   |

---

## `ProteinComplexGroupSplitter`

**Purpose**  
Creates balanced datasets of positive (same-complex) and negative (different-complex) pairs, assigns them to **supergroups** based on network connectivity of complexes, and provides group-aware train/test splitting.

### Parameters
- `positive_label` (default = `1`): Value in `Same_Complex` treated as positive.  
- `negative_ratio` (default = `3`): Ratio of negatives to positives to sample.  
- `random_state` (default = `42`): Random seed for reproducibility.

### Key Methods

#### `_extract_positive_pairs(df)`
- Input: `DataFrame` with `Same_Complex`, `complex_id_A`, `complex_id_B`.  
- Output: Subset of positives plus a `complex_pair` column (`frozenset` of the two complex IDs).

#### `_build_supergroups(pos_df)`
- Builds a graph where complexes are nodes; edges connect complexes that share pairs.  
- Finds connected components = **supergroups**.  
- Assigns each positive pair a `group_id`.  
- Output: `(pos_df_with_group_ids, complex_to_supergroup_dict)`.

#### `_assign_negative_groups(df, pos_group_ids)`
- Samples negatives (`Same_Complex != positive_label`).  
- Size capped at `negative_ratio × (# positives)`.  
- Randomly assigns negatives to one of the positive group IDs.  
- Output: `neg_df` with `group_id`.

#### `prepare(df)`
- Runs full pipeline: extract positives, build supergroups, sample negatives.  
- Returns a DataFrame with both positives and negatives, annotated with `group_id`.

#### `group_split(test_size=0.25)`
- Performs **GroupShuffleSplit**.  
- Ensures train/test splits do not mix group IDs.  
- Output:  
  `(X_train, X_test, y_train, y_test, groups_train, groups_test)`

---

## `stratified_group_split_balanced_per_fold`

**Purpose**  
Performs stratified, group-based cross-validation with balanced classes per fold. Guarantees:

1. Group integrity (`group_col` not split across folds).  
2. Stratification by label.  
3. Equal numbers of positives and negatives in train/test within each fold.

### Parameters
- `df` (`DataFrame`): Dataset with row-level data.  
- `label_col` (`str`): Binary label column (e.g., `Same_Complex_Multi`).  
- `group_col` (`str`): Group identifier (e.g., `pair_id`).  
- `n_splits` (`int`): Number of CV folds.  
- `random_state` (`int`): Random seed.

### Output
- Returns a `list` of `(train_idx, test_idx)` tuples.  
- Each element contains the row indices of train and test sets.  

### Example
```python
folds = stratified_group_split_balanced_per_fold(df, 'Same_Complex_Multi', 'pair_id')
# folds[0] -> (train_indices, test_indices)

