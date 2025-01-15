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


# original functions
def cal_AUROC(sorted_gene_list, positive_set):
    start_time = time.time()
    
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
    
    end_time = time.time()
    print(f"Time (Original): {end_time-start_time:.4f}")
    return AUC

def cal_AUPRC(sorted_gene_list, positive_set):
    start_time = time.time()
    
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
    
    end_time = time.time()
    print(f"Time (Original): {end_time-start_time:.4f}")
    return AUC


# optimize runtime with vectorized computation

def cal_AUROC_vectorized(sorted_gene_list, positive_set):
    """Calculate AUROC with vectorized computation using NumPy arrays.

    Args:
        sorted_gene_list (list): List of identified genes.
        positive_set (set): Set of identified genes in the positive class.

    Returns:
        AUC (float): The AUROC calculated.
    """    
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

def cal_AUPRC_vectorized(sorted_gene_list, positive_set):    
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


def main():
    # print(cpu_count())

    sorted_gene_list, positive_set = create_dummy_inputs(5000, positive_fraction=0.9)
    print(sorted_gene_list)
    print(positive_set)
    
    # AUROC
    print("AUROC:")
    AUROC_1 = cal_AUROC(sorted_gene_list, positive_set)
    print(f"AUROC with original function: {AUROC_1:.4f}")
    start_time = time.time() 
    AUROC_2 = cal_AUROC_vectorized(sorted_gene_list, positive_set)
    end_time = time.time()
    print(f"Time (Vectorized): {end_time-start_time:.4f}")
    print(f"AUROC with vectorized function: {AUROC_2:.4f}")

    # optimize runtime with parallel programming
    # num_pairs = 4
    
    # 1 - multiprocessing
    # processes = []
    
    # for _ in range(pairs): # number loops = number of pairs?
    #     p1 = Process(target=cal_AUROC_vectorized, args=(sorted_gene_list, positive_set))
    #     p2 = Process(target=cal_AUPRC_vectorized)
    #     p1.start()
    #     p2.start()
    #     processes.append(p1)
    #     processes.append(p2)

    # for process in processes:
    #     process.join()

    # 2 - concurrent.futures

    start_time_p = time.time()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(cal_AUROC_vectorized, sorted_gene_list, positive_set)]
        
        for f in concurrent.futures.as_completed(results):
            print(f.result())
    
    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #     results = executor.map(cal_AUROC_vectorized, sorted_gene_list, positive_set)
    
    #     for result in results:
    #         print(result)
    
    end_time_p = time.time()
    
    print(f"Time (Parallel Programming): {end_time_p-start_time_p:.4f}")
    print("----------------------------------------")

    
    # AUPRC
    print("AUPRC:")
    AUPRC_1 = cal_AUPRC(sorted_gene_list, positive_set)
    print(f"AUORC with original function: {AUPRC_1:.4f}")
    AUPRC_2 = cal_AUPRC_vectorized(sorted_gene_list, positive_set)
    print(f"AUPRC with vectorized function: {AUPRC_2:.4f}")
    # AUPRC_3 = cal_AUROC_parallel(sorted_gene_list, positive_set)
    # print(f"AUROC vectorized & run parallel: {AUPRC_3:.4f}")

if __name__ == '__main__':
    main()