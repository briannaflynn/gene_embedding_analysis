#!/usr/bin/env python3
#Date: Dec 05, 2024
#Author: Brianna Flynn
#Description:
# Refactoring calculation.py to optimize for testing multiple thresholds 

import sys
import pickle
import numpy as np
import pandas as pd
import math
from sklearn import metrics
import time
from multiprocessing import Pool, cpu_count
from multiprocessing import shared_memory
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

# main parameters to sweep
# GOBP_term_size_threshold (currently sys.argv[1]), 100 used most recently
# GOBP term level (currently using 2)
# the number of observations necessary for each gene (currently using 1500)

def generate_dummy_data(num_genes=100, num_features=10):
    """
    Generates a dummy DataFrame with random values for testing.
    
    Parameters:
    - num_genes (int): Number of genes (rows).
    - num_features (int): Number of features (columns).
    
    Returns:
    - df (pd.DataFrame): Dummy DataFrame with random values.
    - domain (list): List of random gene names as rows.
    """
    # Create random data
    df = pd.DataFrame(np.random.rand(num_genes, num_features), 
                      index=[f"gene_{i}" for i in range(num_genes)])
    
    # Use all genes as the domain
    domain = df.index.tolist()
    
    return df, domain

def time_function(func, *args, **kwargs):
    """
    Times the execution of a given function.
    
    Parameters:
    - func (callable): The function to time.
    - *args, **kwargs: Arguments to pass to the function.
    
    Returns:
    - result: The output of the function.
    - elapsed_time (float): Execution time in seconds.
    """
    start_time = time.time()
    result = func(*args, **kwargs)
    elapsed_time = time.time() - start_time
    return result, elapsed_time

def compute_pairwise_spearman_loop(df, domain):
    """
    Computes pairwise Spearman correlations using a nested loop.
    
    Parameters:
    - df (pd.DataFrame): Input DataFrame.
    - domain (list): List of genes to consider.
    
    Returns:
    - dict: Dictionary of gene pairs and their Spearman correlations.
    """
    gene_gene2_corr = {}
    for gene in domain:
        for gene2 in domain:
            if gene2 == gene:
                continue
            gene_, gene2_ = sorted([gene, gene2])
            if (gene_, gene2_) in gene_gene2_corr:
                continue
            cor_value = df.loc[gene_].corr(df.loc[gene2_], method="spearman")
            gene_gene2_corr[(gene_, gene2_)] = cor_value
    return gene_gene2_corr

def compute_correlations_for_chunk(chunk, filtered_df, row_indices, k):
    """
    Computes correlations for a subset of row pairs.
    
    Parameters:
    - chunk (list of tuple): List of row index pairs (i, j).
    - filtered_df (pd.DataFrame): Filtered DataFrame with relevant rows.
    - row_indices (list): List of row names in the filtered DataFrame.
    - k (int): Minimum offset for the upper triangle (default k=1).
    
    Returns:
    - list of tuple: Computed correlations in the form (gene1, gene2, correlation).
    """
    result = []
    for i, j in chunk:
        gene1 = row_indices[i]
        gene2 = row_indices[j]
        correlation = filtered_df.loc[gene1].corr(filtered_df.loc[gene2], method='spearman')
        result.append((gene1, gene2, correlation))
    return result

def compute_pairwise_spearman_parallel(df, domain, k=1, n_jobs=None):
    """
    Computes pairwise Spearman correlations for the specified rows of a DataFrame
    in parallel and returns the upper triangle of the correlation matrix as a long-form DataFrame.

    Parameters:
    - df (pd.DataFrame): The input DataFrame (genes as rows, features as columns).
    - domain (list or set): The list of row indices to filter for pairwise correlation.
    - k (int): Minimum offset for the upper triangle (default k=1 excludes the diagonal).
    - n_jobs (int): Number of parallel jobs (default None uses all available CPUs).

    Returns:
    - pd.DataFrame: A long-form DataFrame with the pairwise correlations.
    """
    # Filter the DataFrame for the specified domain
    filtered_df = df.loc[domain]
    row_indices = list(filtered_df.index)
    num_rows = len(row_indices)

    # Generate pairs for the upper triangle
    all_pairs = [(i, j) for i in range(num_rows) for j in range(i + k, num_rows)]

    # Divide the pairs into chunks for parallel processing
    n_jobs = n_jobs or cpu_count()  # Use all available CPUs if n_jobs is not specified
    chunk_size = (len(all_pairs) + n_jobs - 1) // n_jobs
    chunks = [all_pairs[i:i + chunk_size] for i in range(0, len(all_pairs), chunk_size)]
    #print(len(chunks))
    print('beginning parallel process')
    # Use multiprocessing to compute correlations in parallel
    with Pool(processes=n_jobs) as pool:
        results = pool.starmap(
            compute_correlations_for_chunk, 
            [(chunk, filtered_df, row_indices, k) for chunk in chunks]
        )

    # Flatten results
    flattened_results = [item for sublist in results for item in sublist]
    #print(flattened_results)

    # Convert to a DataFrame
    pairwise_correlations = pd.DataFrame(flattened_results, columns=['Row', 'Column', 'Correlation'])
    #print(pairwise_correlations)

    return pairwise_correlations

def compute_correlations_for_chunk_shared(chunk, shm_name, shape, row_indices):
    """
    Computes Spearman correlations for a subset of row pairs using shared memory.

    Parameters:
    - chunk (list of tuple): List of row index pairs (i, j).
    - shm_name (str): Name of the shared memory block.
    - shape (tuple): Shape of the shared data array.
    - row_indices (list): List of row indices corresponding to the dataset.

    Returns:
    - list of tuple: Computed correlations in the form (gene1, gene2, correlation).
    """
    # Attach to the shared memory block
    shm = shared_memory.SharedMemory(name=shm_name)
    shared_data = np.ndarray(shape, buffer=shm.buf)
    
    results = []
    for i, j in chunk:
        gene1 = row_indices[i]
        gene2 = row_indices[j]
        # Compute Spearman correlation
        correlation, _ = spearmanr(shared_data[i], shared_data[j])
        results.append((gene1, gene2, correlation))
    
    # Close the shared memory block in the worker
    shm.close()
    return results

def compute_pairwise_spearman_parallel_shared(df, domain, k=1, n_jobs=None):
    """
    Computes pairwise Spearman correlations for the specified rows of a DataFrame
    in parallel using shared memory, returning the upper triangle as a long-form DataFrame.

    Parameters:
    - df (pd.DataFrame): The input DataFrame (genes as rows, features as columns).
    - domain (list or set): The list of row indices to filter for pairwise correlation.
    - k (int): Minimum offset for the upper triangle (default k=1 excludes the diagonal).
    - n_jobs (int): Number of parallel jobs (default None uses all available CPUs).

    Returns:
    - pd.DataFrame: A long-form DataFrame with the pairwise correlations.
    """
    # Filter the DataFrame for the specified domain
    filtered_df = df.loc[domain]
    row_indices = list(filtered_df.index)
    num_rows = len(row_indices)

    # Create shared memory
    shm = shared_memory.SharedMemory(create=True, size=filtered_df.values.nbytes)
    shared_array = np.ndarray(filtered_df.shape, dtype=filtered_df.values.dtype, buffer=shm.buf)
    shared_array[:] = filtered_df.values  # Copy data to shared memory

    # Generate pairs for the upper triangle
    all_pairs = [(i, j) for i in range(num_rows) for j in range(i + k, num_rows)]

    # Divide the pairs into chunks for parallel processing
    #n_jobs = n_jobs or len(all_pairs)  # fix this - doesn't get cpu count
    n_jobs = n_jobs or cpu_count()
    chunk_size = (len(all_pairs) + n_jobs - 1) // n_jobs
    chunks = [all_pairs[i:i + chunk_size] for i in range(0, len(all_pairs), chunk_size)]

    try:
        # Use multiprocessing to compute correlations in parallel
        with Pool(processes=n_jobs) as pool:
            results = pool.starmap(
                compute_correlations_for_chunk_shared, 
                [(chunk, shm.name, filtered_df.shape, row_indices) for chunk in chunks]
            )
    finally:
        # Cleanup shared memory
        shm.close()
        shm.unlink()

    # Flatten results
    flattened_results = [item for sublist in results for item in sublist]

    # Convert to a DataFrame
    pairwise_correlations = pd.DataFrame(flattened_results, columns=["Row", "Column", "Correlation"])

    return pairwise_correlations


def benchmark_runtime(func, df, domain, fractions):
    """
    Benchmarks runtime of a function at different scales and fits a predictive model.
    
    Parameters:
    - func (callable): The function to benchmark.
    - df (pd.DataFrame): The full input DataFrame.
    - domain (list): Full list of genes (rows) to consider.
    - fractions (list): Fractions of the full data to use for benchmarking (e.g., [0.1, 0.25, 0.5]).
    
    Returns:
    - times (list): Measured runtimes for each fraction.
    - sizes (list): Sizes of the domains for each fraction.
    """
    times = []
    sizes = []
    
    for fraction in fractions:
        print('starting fraction', fraction)
        sample_size = int(len(domain) * fraction)
        sampled_domain = domain[:sample_size]
        
        start_time = time.time()
        func(df, sampled_domain)
        elapsed_time = time.time() - start_time
        
        times.append(elapsed_time)
        sizes.append(sample_size)
    
    return sizes, times

def fit_and_predict(sizes, times, target_sizes, degree=2):
    """
    Fits a polynomial regression model and predicts runtimes for target sizes.
    
    Parameters:
    - sizes (list): Sizes of the domains used for benchmarking.
    - times (list): Measured runtimes.
    - target_sizes (list): Sizes to predict runtimes for.
    - degree (int): Degree of the polynomial to fit.
    
    Returns:
    - predictions (list): Predicted runtimes for the target sizes.
    """
    # Reshape data for regression
    sizes = np.array(sizes).reshape(-1, 1)
    times = np.array(times).reshape(-1, 1)
    target_sizes = np.array(target_sizes).reshape(-1, 1)
    
    # Fit polynomial regression model
    poly = PolynomialFeatures(degree=degree)
    sizes_poly = poly.fit_transform(sizes)
    model = LinearRegression()
    model.fit(sizes_poly, times)
    
    # Predict runtimes for target sizes
    target_sizes_poly = poly.transform(target_sizes)
    predictions = model.predict(target_sizes_poly)
    
    return predictions.flatten()


# EXECUTION
if __name__ == "__main__":
    df = pd.read_csv("ratio_of_n_and_n_cells_cell_type.gene.tissue-celltype.ver2.tsv", sep="\t", header=0, index_col=0)
    domain = pd.read_csv("all_proteins_considered.intersection_df_human-ref-proteome.csv")
    domain = domain['Domain'].to_list()
    print(df.loc[domain].shape)
    
    # Compute pairwise Spearman correlations using shared memory
    parallel_results = compute_pairwise_spearman_parallel_shared(df, domain, k=1, n_jobs=4)
    
    print(parallel_results.head())

    #### PERFORMANCE TEST ####
    # df, domain = generate_dummy_data(num_genes=20000, num_features=3500)  # Generate dummy data
    # fractions = [0.001, 0.005, 0.01]  # Benchmark at 10%, 25%, 50% scales
    # sizes, times = benchmark_runtime(compute_pairwise_spearman_loop, df, domain, fractions)
    
    # # Predict runtimes for larger sizes
    # target_sizes = [12000, 15000, 20000]  # Sizes to project
    # predicted_times = fit_and_predict(sizes, times, target_sizes)


