import re
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import (
    roc_curve, auc,
    precision_recall_curve, average_precision_score
)

"""
Example
-------
INPUT - final_df, xf, X_full

X_full = pd.read_pickle("/home/ubuntu/full_dataset.pkl")
with open("model_results.pkl", "rb") as f:
    _all_model_results_ = pickle.load(f)

for mod in list(_all_model_results_.keys()):
    for p in ['predictions', 'probabilities', 'decision_scores']:
        col_name = f'{mod}_{p}'
        X_full[col_name] = _all_model_results_[mod][p]

xf = final_df[final_df['Test'] == True][['GeneAB']].merge(X_full, on="GeneAB")

plot_roc_curves(xf, "Same_Complex")
plot_pr_curves(xf, "Same_Complex")
"""


def get_model_prefixes(df):
    """
    Extract unique model prefixes from DataFrame columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing model output columns such as
        <model>_probabilities or <model>_decision_scores.

    Returns
    -------
    list of str
        Sorted list of model prefixes.
    """
    prefixes = set()
    for col in df.columns:
        match = re.match(r"(.+)(_probabilities|_decision_scores)$", col)
        if match:
            prefixes.add(match.group(1))
    return sorted(prefixes)


def get_valid_scores(df, prefix):
    """
    Retrieve a valid score array for a given model prefix.

    Prefers probabilities if available, otherwise falls back to decision scores.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing model outputs.
    prefix : str
        Model prefix string.

    Returns
    -------
    np.ndarray or None
        Array of scores if available, otherwise None.
    """
    for suffix in ["_probabilities", "_decision_scores"]:
        colname = f"{prefix}{suffix}"
        if colname not in df.columns:
            continue

        y_score = df[colname].values
        if all(v is None for v in y_score):
            continue
        try:
            return np.array(y_score, dtype=float)
        except Exception:
            continue

    return None


def plot_roc_curves(df, y_col, savepath="full_train_test_roc.png"):
    """
    Plot ROC curves for all models with valid scores.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing labels and model outputs.
    y_col : str
        Column name of ground truth labels.
    savepath : str, default="full_train_test_roc.png"
        Path to save the ROC curve plot.
    """
    y_true = df[y_col].values
    prefixes = get_model_prefixes(df)

    plt.figure(figsize=(10, 7))
    for prefix in prefixes:
        y_score = get_valid_scores(df, prefix)
        if y_score is None:
            print(f"Skipping {prefix} (no valid scores)")
            continue

        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, label=f"{prefix} (AUC = {roc_auc:.2f})")

    plt.plot([0, 1], [0, 1], "k--", lw=1)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curve")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(savepath, dpi=300)
    plt.show()


def plot_pr_curves(df, y_col, savepath="full_train_test_auprc.png"):
    """
    Plot Precision-Recall curves for all models with valid scores.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing labels and model outputs.
    y_col : str
        Column name of ground truth labels.
    savepath : str, default="full_train_test_auprc.png"
        Path to save the PR curve plot.
    """
    y_true = df[y_col].values
    prefixes = get_model_prefixes(df)

    plt.figure(figsize=(10, 7))
    for prefix in prefixes:
        y_score = get_valid_scores(df, prefix)
        if y_score is None:
            print(f"Skipping {prefix} (no valid scores)")
            continue

        precision, recall, _ = precision_recall_curve(y_true, y_score)
        ap = average_precision_score(y_true, y_score)
        plt.plot(recall, precision, label=f"{prefix} (AP = {ap:.2f})")

    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("Precision-Recall Curve")
    plt.legend(loc="lower left")
    plt.tight_layout()
    plt.savefig(savepath, dpi=300)
    plt.show()
