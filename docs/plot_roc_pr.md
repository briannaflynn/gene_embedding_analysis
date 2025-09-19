# Model Evaluation Utilities

`ROC_PR.py`

This module provides helper functions for extracting model scores and visualizing performance across multiple classifiers. It supports plotting **ROC curves** and **Precision–Recall curves** using predictions stored in a `pandas.DataFrame`.

---

## Example Workflow

```python
X_full = pd.read_pickle("/home/ubuntu/full_dataset.pkl")
with open("model_results.pkl", "rb") as f:
    _all_model_results_ = pickle.load(f)

# Add model outputs to the feature DataFrame
for mod in list(_all_model_results_.keys()):
    for p in ['predictions', 'probabilities', 'decision_scores']:
        col_name = f'{mod}_{p}'
        X_full[col_name] = _all_model_results_[mod][p]

# Select only test rows and join with model outputs
xf = final_df[final_df['Test'] == True][['GeneAB']].merge(X_full, on="GeneAB")

# Plot metrics
plot_roc_curves(xf, "Same_Complex")
plot_pr_curves(xf, "Same_Complex")
```
---

## `get_model_prefixes(df)`

### Purpose
Extracts the set of unique model prefixes (e.g., "RandomForest", "LogisticRegression") from the DataFrame column names.

### Input
- `df` (`pd.DataFrame`): Must contain model output columns of the form `<model>_probabilities` or `<model>_decision_scores`.

### Output
- `list[str]`: Sorted list of unique model prefixes.

---

## `get_valid_scores(df, prefix)`

### Purpose
Fetches a valid score array for a given model prefix:
- Prefers `_probabilities` if present.
- Falls back to `_decision_scores` if probabilities are unavailable.

### Input
- `df` (`pd.DataFrame`): DataFrame containing model outputs.
- `prefix` (`str`): Model name prefix.

### Output
- `np.ndarray` of floats if valid scores are found.
- `None` if neither probabilities nor decision scores are available.

---

## `plot_roc_curves(df, y_col, savepath="full_train_test_roc.png")`

### Purpose
Generates and saves a **Receiver Operating Characteristic (ROC)** curve for each model with valid scores.

### Input
- `df` (`pd.DataFrame`): DataFrame containing ground truth labels and model outputs.
- `y_col` (`str`): Column name with binary ground truth labels (e.g., "Same_Complex").
- `savepath` (`str`, default = "full_train_test_roc.png"): Path to save the plot.

### Output
- A saved PNG file containing the ROC curves.
- A matplotlib figure showing:
  - ROC curves for each model.
  - Diagonal line for random baseline.
  - AUC values in the legend.

---

## `plot_pr_curves(df, y_col, savepath="full_train_test_auprc.png")`

### Purpose
Generates and saves a **Precision–Recall (PR)** curve for each model with valid scores.

### Input
- `df` (`pd.DataFrame`): DataFrame containing ground truth labels and model outputs.
- `y_col` (`str`): Column name with binary ground truth labels.
- `savepath` (`str`, default = "full_train_test_auprc.png"): Path to save the plot.

### Output
- A saved PNG file containing the PR curves.
- A matplotlib figure showing:
  - Precision vs. Recall for each model.
  - Average Precision (AP) values in the legend.

---

## Expected DataFrame Schema

Your input `df` must contain:

- **Ground truth column**: e.g., `Same_Complex` (binary labels).  
- **Model output columns**: For each model `<prefix>`, at least one of:
  - `<prefix>_probabilities` (preferred).  
  - `<prefix>_decision_scores` (fallback).

Optional:
- `<prefix>_predictions` can be included but is not used in plotting.

---

## Example Output Plots

- **ROC Curve**: Shows trade-off between sensitivity (TPR) and false positive rate (FPR).  
- **PR Curve**: Shows trade-off between precision and recall, especially informative in imbalanced datasets.