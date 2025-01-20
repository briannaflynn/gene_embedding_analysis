import matplotlib.pyplot as plt
import pandas as pd

combined_df = pd.read_csv('Combined_Data_for_Top1AUPR_and_Protein_Coverage.csv')
sorted_df = combined_df.sort_values(by='Number of Proteins')

x_sorted = sorted_df['Number of Proteins']
y_sorted = sorted_df['Top1AUPRC']
y_secondary = sorted_df['Top1AUROC']

def double_axis_1():
    fig, ax1 = plt.subplots(figsize=(12, 6))
    ax1.plot(x_sorted, y_sorted, color='darkblue', label='Top1AUPR', marker='o', linestyle='--')
    ax1.set_xlabel("Number of Proteins")
    ax1.set_ylabel("Top1AUPR", color='black')
    ax1.tick_params(axis='y', labelcolor='black')
    ax1.grid(True, linestyle='--', alpha=0.7)

    ax2 = ax1.twinx()
    ax2.plot(x_sorted, y_secondary, color='orange', label='Top1AUROC', marker='o', linestyle='--')
    ax2.set_ylabel("Top1AUROC", color='black')
    ax2.tick_params(axis='y', labelcolor='black')
    ax2.set_ylim(ax1.get_ylim())

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    fig.legend(lines1 + lines2, labels1 + labels2, loc='upper right', bbox_to_anchor=(0.9, 0.9))

    fig.suptitle("Top1AUPR and Top1AUROC vs Number of Proteins Remaining after Threshold Filter", fontsize=14)
    plt.show()

def double_axis_2():
    fig, ax1 = plt.subplots(figsize=(12, 6))
    ax1.plot(x_sorted, y_sorted, color='darkblue', label='Top1AUPR', marker='o', linestyle='--')
    ax1.set_xlabel("Number of Proteins")
    ax1.set_ylabel("Top1AUPR", color='black')
    ax1.tick_params(axis='y', labelcolor='black')
    ax1.grid(True, linestyle='--', alpha=0.7)

    ax2 = ax1.twinx()
    ax2.plot(x_sorted, y_secondary, color='orange', label='Top1AUROC', marker='o', linestyle='--')
    ax2.set_ylabel("Top1AUROC", color='black')
    ax2.tick_params(axis='y', labelcolor='black')

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    fig.legend(lines1 + lines2, labels1 + labels2, loc='upper right', bbox_to_anchor=(0.9, 0.9))

    fig.suptitle("Top1AUPR and Top1AUROC vs Number of Proteins Remaining after Threshold Filter", fontsize=14)
    plt.show()

def single_line_plot():
    plt.figure(figsize=(12, 6))
    plt.plot(x_sorted, y_sorted, color='darkblue', marker='o', linestyle='--', label='Top1AUPR vs Number of Proteins')
    plt.title("Top1AUPR vs Number of Proteins as a Function of Threshold")
    plt.xlabel("Number of Proteins")
    plt.ylabel("Top1AUPR")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.show()

def double_axis_3():
    fig, ax1 = plt.subplots(figsize=(12, 6))
    ax1.plot(sorted_df['Threshold (Observations)'], sorted_df['Number of Proteins'], color='darkgreen', label='Number of Proteins', marker='o', linestyle='--')
    ax1.set_xlabel("Threshold (Observations)")
    ax1.set_ylabel("Number of Proteins", color='darkgreen')
    ax1.tick_params(axis='y', labelcolor='darkgreen')
    ax1.grid(True, linestyle='--', alpha=0.7)

    ax2 = ax1.twinx()
    ax2.plot(sorted_df['Threshold (Observations)'], sorted_df['Top1AUPRC'], color='darkblue', label='Top1AUPR', marker='o', linestyle='--')
    ax2.set_ylabel("Top1AUPR", color='darkblue')
    ax2.tick_params(axis='y', labelcolor='darkblue')

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    fig.legend(lines1 + lines2, labels1 + labels2, loc='upper right', bbox_to_anchor=(0.9, 0.9))

    fig.suptitle("Top1AUPR and Number of Proteins Remaining vs Threshold (Observations)", fontsize=14)
    plt.show()
