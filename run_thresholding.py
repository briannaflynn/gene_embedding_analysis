from merge_pairwise_with_gobp import *
from calculate_benchmark_metrics import *
import pickle 
import sys

DIR = "/scratch/05515/bflynn/gene_embedding_benchmark/"
print("Reading in pairwise")
pairdf =  pd.read_pickle(DIR + 'pairwise_correlations.pkl')
df_rename_by_idx(pairdf)

print("Reading in gobp")
gobpdf =  pd.read_pickle(DIR + 'gobp_size100_level2_debug.pkl')
df_rename_by_idx(gobpdf)

def save_dict_to_pickle(dictionary, filename):

    try:
        with open(filename, 'wb') as file:
            pickle.dump(dictionary, file)
        print(f"Dictionary successfully saved to {filename}")
    except Exception as e:
        print(f"An error occurred while saving the dictionary: {e}")

def runner(pairwise, gobp, filter_num):

    print(f'Starting analysis for: all_proteins_considered.intersection_df_human-ref-proteome_{filter_num}.csv')

    filter_obs_df = pd.read_csv(DIR + f'all_proteins_considered.intersection_df_human-ref-proteome_{filter_num}.csv')
    print(filter_obs_df.shape)
    df_pair = filter_by_obs(pairwise, filter_obs_df)
    print(df_pair.shape)
    print('Startng merge')

    merged = merge_pairwise_dataframes(df_pair, gobp) 
    print(merged.shape)

    calculator = MetricCalculatorSKLearn()

    print('Calculate metrics including top-k (default 1%) AUPRC and AUROC')
    results = calculator.calculate_metrics(merged, score_col="Correlation", target_col="GOBP_value")
    
    # Print Results
    print(f"Results for all_proteins_considered.intersection_df_human-ref-proteome_{filter_num}")
    print(f"AUROC: {results['AUROC']}")
    print(f"AUPRC: {results['AUPRC']}")
    print(f"Top-k (1%) AUPRC: {results['TopKAUPRC']}")
    print(f"Top-k (1%) AUROC: {results['TopKAUROC']}")

    save_dict_to_pickle(results, DIR + f'pairwise_corr.gobp100.all_proteins_considered.intersection_df_human-ref-proteome_{filter_num}.pkl')

if __name__ == "__main__":
    
    start_number = 150  # Replace with your desired starting number
    end_number = 1500    # Replace with your desired ending number

    # Loop to increase the number by 50 each time
    current_number = start_number
    while current_number <= end_number:
        runner(pairdf, gobpdf, current_number)
        current_number += 50