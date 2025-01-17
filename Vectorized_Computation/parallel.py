import numpy as np
from sklearn import metrics
import random
import time
from multiprocessing import cpu_count, Pool
import concurrent.futures

def create_dummy_inputs(list_length, positive_fraction=0.1, seed=42):
    """
    Create dummy inputs for testing cal_AUROC and cal_AUPRC.
    
    Parameters:
    - list_length (int): Length of the sorted_gene_list.
    - positive_fraction (float): Fraction of the genes to include in the positive set.
    - seed (int): Random seed for reproducibility.
    
    Returns:
    - sorted_gene_list (list): List of dummy gene identifiers.
    - positive_set (set): Set of positive dummy gene identifiers.
    """
    random.seed(seed)
    sorted_gene_list = [f"gene_{i}" for i in range(list_length)]
    num_positive = int(list_length * positive_fraction)
    positive_set = set(random.sample(sorted_gene_list, num_positive))
    
    return sorted_gene_list, positive_set

# Example usage
# sorted_gene_list, positive_set = create_dummy_inputs(list_length=10000, positive_fraction=0.2)


# 1 - original functions

def cal_AUROC(sorted_gene_list, positive_set):
    total = len(sorted_gene_list)
    actual_positive = positive_set & set(sorted_gene_list)
    # x: FPR = FP / Actual Negative
    # y: TPR = TP / Actual Positive
    AP = len(actual_positive)
    AN = total - AP
    FP = 0
    TP = 0
    coordinates = []
    for gene in sorted_gene_list:
        if gene in actual_positive:
            TP += 1
        else:
            FP += 1
        coordinates.append((FP/AN, TP/AP))
    FPR = [x[0] for x in coordinates]
    TPR = [x[1] for x in coordinates]   
    AUC = metrics.auc(FPR, TPR)
    return AUC

def cal_AUPRC(sorted_gene_list, positive_set):    
    total = len(sorted_gene_list)
    actual_positive = positive_set & set(sorted_gene_list)
    # x: recall = (TP / Actual Positive)
    # y: precision = (TP / Predicted Positive)
    AP = len(actual_positive)
    PP = 0
    TP = 0
    coordinates = []
    for gene in sorted_gene_list:
        PP += 1
        if gene in actual_positive:
            TP += 1
        coordinates.append((TP/AP, TP/PP))
    recall = [x[0] for x in coordinates]
    precision = [x[1] for x in coordinates]
    AUC = metrics.auc(recall, precision)
    return AUC


# 2 - optimize runtime with vectorized computation

def cal_AUROC_vectorized(input):
    sorted_gene_list, positive_set = input
    total = len(sorted_gene_list)
    actual_positive = positive_set & set(sorted_gene_list)
    # x: FPR = FP / Actual Negative
    # y: TPR = TP / Actual Positive
    AP = len(actual_positive)
    AN = total - AP
    positive_array = np.array([1 if gene in actual_positive else 0 for gene in sorted_gene_list])
    TP_array = np.cumsum(positive_array)
    FP_array = np.cumsum(1-positive_array)
    FPR = FP_array / AN
    TPR = TP_array / AP
    AUC = metrics.auc(FPR, TPR)
    return AUC

def cal_AUPRC_vectorized(input):
    sorted_gene_list, positive_set = input
    total = len(sorted_gene_list)
    actual_positive = positive_set & set(sorted_gene_list)
    # x: recall = (TP / Actual Positive)
    # y: precision = (TP / Predicted Positive)
    AP = len(actual_positive)
    AN = total - AP
    PP = len(sorted_gene_list)
    positive_array = np.array([1 if gene in actual_positive else 0 for gene in sorted_gene_list])
    TP_array = np.cumsum(positive_array)
    PP_array = np.arange(1, PP+1, 1)
    recall = TP_array / AP
    precision = TP_array / PP_array
    AUC = metrics.auc(recall, precision)
    return AUC

# 3 - optimize runtime with parallel processing

def cal_AUROC_parallel(inputs, n_jobs=None):
    n_jobs = n_jobs or cpu_count()

    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = executor.map(cal_AUROC_vectorized, inputs)
    
        # for result in results:
        #     print(result)
    return results

def cal_AUPRC_parallel(inputs, n_jobs=None):
    n_jobs = n_jobs or cpu_count()

    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = executor.map(cal_AUPRC_vectorized, inputs)

        # for result in results:
        #     print(result)
    return results


def main():
    # print(cpu_count())

    sorted_gene_list, positive_set = create_dummy_inputs(100000, positive_fraction=0.9)
    # print(sorted_gene_list)
    # print(positive_set)

    num_pairs = 40
    
    # AUROC
    print("AUROC:")

    # 1 - original
    start_time_roc_1 = time.time() 
    AUROC_1 = cal_AUROC(sorted_gene_list, positive_set)
    end_time_roc_1 = time.time()
    time_original_AUROC = (end_time_roc_1-start_time_roc_1) * num_pairs
    print(f"Time (Original): {time_original_AUROC:.4f}")
    print(f"AUROC with original function: {AUROC_1:.4f}")

    # 2 - vectorized
    start_time_roc_2 = time.time() 
    AUROC_2 = cal_AUROC_vectorized((sorted_gene_list, positive_set))
    end_time_roc_2 = time.time()
    print(f"Time (Vectorized): {(end_time_roc_2-start_time_roc_2) * num_pairs:.4f}")
    print(f"AUROC with vectorized function: {AUROC_2:.4f}")

    # 3 - parallel programming
    inputs = []
    for i in range(num_pairs):
        sorted_list, positive_set = create_dummy_inputs(5000, positive_fraction=0.9, seed=42+i)
        inputs.append((sorted_list, positive_set))

    start_time_roc_3 = time.time()
    results_AUROC = cal_AUROC_parallel(inputs)
    end_time_roc_3 = time.time()
    time_vectorized_parallel_AUROC = end_time_roc_3-start_time_roc_3
    print(f"Time (Parallel Programming): {time_vectorized_parallel_AUROC:.4f}")
    l_roc = list(results_AUROC)
    avg_result_roc = sum(l_roc) / len(l_roc)
    print(f"AUROC with parallel programming on vectorized function: {avg_result_roc:.4f}")

    print("----------------------------------------")
    
    # AUPRC
    print("AUPRC:")
    
    # 1 - original
    start_time_prc_1 = time.time() 
    AUPRC_1 = cal_AUPRC(sorted_gene_list, positive_set)
    end_time_prc_1 = time.time()
    time_original_AUPRC = (end_time_prc_1-start_time_prc_1) * num_pairs
    print(f"Time (Original): {time_original_AUPRC:.4f}")
    print(f"AUPRC with original function: {AUPRC_1:.4f}")

    # 2 - vectorized
    start_time_prc_2 = time.time() 
    AUPRC_2 = cal_AUPRC_vectorized((sorted_gene_list, positive_set))
    end_time_prc_2 = time.time()
    print(f"Time (Vectorized): {(end_time_prc_2-start_time_prc_2) * num_pairs:.4f}")
    
    print(f"AUPRC with vectorized function: {AUPRC_2:.4f}")

    # 3 - parallel programming
    start_time_prc_3 = time.time()
    results_AUPRC = cal_AUPRC_parallel(inputs)
    end_time_prc_3 = time.time()
    time_vectorized_parallel_AUPRC = end_time_prc_3-start_time_prc_3
    print(f"Time (Parallel Programming): {time_vectorized_parallel_AUPRC:.4f}")
    l_prc = list(results_AUPRC)
    avg_result_prc = sum(l_prc) / len(l_prc)
    print(f"AUPRC with parallel programming on vectorized function: {avg_result_prc:.4f}")

    
    # compare runtime
    speedup_AUROC = (time_original_AUROC * num_pairs) / time_vectorized_parallel_AUROC
    print(f'- AUROC Speedup (Vectorized / Parallel vs. Original): {speedup_AUROC}x')

    speedup_AUPRC = (time_original_AUPRC * num_pairs) / time_vectorized_parallel_AUPRC
    print(f'- AUPRC Speedup (Vectorized / Parallel vs. Original): {speedup_AUPRC}x')


if __name__ == '__main__':
    main()