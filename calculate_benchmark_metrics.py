#!/usr/bin/env python3
#Date: Dec 09, 2024
#Author: Brianna Flynn
#Description:
#This script takes data formatted from merge_pairwise_with_gobp.py and gets AUROC, AUPRC, and top 1% AUROC and AUPRC metrics as well as vectors for curve plotting - export as JSON not csv (see last line), probably should be pkl since files are large
import pandas as pd
from sklearn.metrics import roc_curve, precision_recall_curve, roc_auc_score, auc
import numpy as np
import json

class MetricCalculatorSKLearn:
    def __init__(self):
        self.results = {}

    def calculate_metrics(self, dataframe, score_col, target_col):
        """
        Calculate AUROC, AUPRC, and top-k metrics using scikit-learn and store intermediate values.
        """
        scores = dataframe[score_col].values
        targets = dataframe[target_col].values

        # Calculate AUROC and ROC curve
        fpr, tpr, _ = roc_curve(targets, scores)
        auroc = roc_auc_score(targets, scores)

        # Calculate Precision-Recall Curve and AUPRC
        precision, recall, _ = precision_recall_curve(targets, scores)
        auprc = auc(recall, precision)

        # Calculate top-k metrics (1% by default)
        top_k_results = self._calculate_top_k(scores, targets, k=0.01)

        # Save results
        self.results = {
            'AUROC': auroc,
            'AUPRC': auprc,
            'TopKAUPRC': top_k_results['AUPRC'],
            'TopKAUROC': top_k_results['AUROC'],
            'FPR': fpr.tolist(),
            'TPR': tpr.tolist(),
            'Precision': precision.tolist(),
            'Recall': recall.tolist()
        }

        return self.results

    def _calculate_top_k(self, scores, targets, k=0.01):
        """
        Calculate the AUPRC and AUROC for the top-k% of scores.

        Args:
        - scores: Array of scores (e.g., cosine similarity or correlation).
        - targets: Array of binary target values (0 or 1).
        - k: Proportion of the top scores to consider (default: 0.01 for 1%).

        Returns:
        - A dictionary with the AUPRC and AUROC for the top-k% scores.
        """
        # Get the threshold for the top-k% scores
        top_k_threshold = np.percentile(scores, 100 * (1 - k))
        top_k_indices = scores >= top_k_threshold

        # Filter scores and targets to the top-k%
        filtered_scores = scores[top_k_indices]
        filtered_targets = targets[top_k_indices]

        if len(filtered_targets) == 0 or sum(filtered_targets) == 0:
            return {'AUPRC': np.nan, 'AUROC': np.nan}  # Handle edge cases

        # Calculate Precision-Recall Curve and AUPRC for the filtered data
        precision, recall, _ = precision_recall_curve(filtered_targets, filtered_scores)
        auprc = auc(recall, precision)

        # Calculate ROC curve and AUROC for the filtered data
        fpr, tpr, _ = roc_curve(filtered_targets, filtered_scores)
        auroc = roc_auc_score(filtered_targets, filtered_scores)

        return {'AUPRC': auprc, 'AUROC': auroc}

    def export_to_csv(self, filepath):
        """
        Export detailed intermediate results for AUROC (FPR, TPR), AUPRC (Precision, Recall),
        and top-k metrics to a CSV file.

        Args:
        - filepath: Path to the CSV file to save results.
        """
        # Gather data for the CSV
        data = {
            'FPR': self.results['FPR'],
            'TPR': self.results['TPR'],
            'Precision': self.results['Precision'],
            'Recall': self.results['Recall']
        }

        # Convert to DataFrame
        df = pd.DataFrame(data)

        # Add overall metrics (AUROC, AUPRC, top-k AUPRC, top-k AUROC) as separate columns
        df['AUROC'] = self.results['AUROC']
        df['AUPRC'] = self.results['AUPRC']

        df['TopKAUPRC'] = self.results['TopKAUPRC']
        df['TopKAUROC'] = self.results['TopKAUROC']

        # Save to CSV
        df.to_csv(filepath, index=False)
        print(f"Detailed results exported to {filepath}")

# data = {
#     "gene1": ["gene_0", "gene_1", "gene_2", "gene_3", "gene_4", "gene_5", "gene_6", "gene_7", "gene_8", "gene_9"],
#     "gene2": ["gene_1", "gene_2", "gene_3", "gene_4", "gene_5", "gene_6", "gene_7", "gene_8", "gene_9", "gene_10"],
#     "correlation": [0.374540, 0.950714, 0.731994, 0.598658, 0.156019, 0.155995, 0.058084, 0.866176, 0.601115, 0.708073],
#     "gobp_value": [0, 0, 1, 0, 1, 1, 0, 1, 0, 0]
# }
# import sys
# # Create the DataFrame


# pickledf = pd.read_pickle('pairwise_150_gobp_100_lvl2_debug.pkl')
# print(pickledf.columns)

# calculator = MetricCalculatorSKLearn()

# # Calculate metrics including top-k (default 1%) AUPRC and AUROC
# results = calculator.calculate_metrics(pickledf, score_col="Correlation", target_col="GOBP_value")

# # Print Results
# print(
# print(f"AUROC: {results['AUROC']}")
# print(f"AUPRC: {results['AUPRC']}")
# print(f"Top-k (1%) AUPRC: {results['TopKAUPRC']}")
# print(f"Top-k (1%) AUROC: {results['TopKAUROC']}")

# print(results)
# # Export results to CSV
# with open('pairwise_150_gobp_100_lvl2.metrics_detailed_results.json', 'w') as fp:
#     json.dump(results, fp)


