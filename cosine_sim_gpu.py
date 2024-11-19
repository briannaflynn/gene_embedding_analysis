import torch
from torch.nn.functional import normalize
import pandas as pd
from itertools import product
import pickle
import argparse
from tqdm import tqdm
import numpy as np
from itertools import islice
import time
from multiprocessing import Pool
import os

def ensure_tensor(embedding_dict):
    """Ensure all embeddings in the dictionary are PyTorch tensors."""
    for key, value in embedding_dict.items():
        if isinstance(value, np.ndarray):
            embedding_dict[key] = torch.tensor(value, dtype=torch.float32)
        elif not isinstance(value, torch.Tensor):
            raise TypeError(f"Embedding for key '{key}' must be either a NumPy array or PyTorch tensor.")
    return embedding_dict


def load_embeddings_from_pkl(pkl_file):
    with open(pkl_file, 'rb') as f:
        embedding_dict = pickle.load(f)
    return ensure_tensor(embedding_dict)

"""Serial write out"""
def write_similarity_table_incrementally(output_file, label_pairs, similarity_values, chunk_size=1000):
    # write similarity table incrementally to the output file since pandas taking a while
    assert len(label_pairs) == len(similarity_values), "Mismatch in label pairs and similarity values lengths."
    
    with open(output_file, 'w') as f:
        # header
        f.write("gene_1,gene_2,cosine_similarity\n")

        # Initialize the tqdm progress bar
        with tqdm(total=len(label_pairs), desc="Writing similarity table", unit="pairs") as pbar:
            for start_idx in range(0, len(label_pairs), chunk_size):
                end_idx = min(start_idx + chunk_size, len(label_pairs))
                chunk_pairs = islice(label_pairs, start_idx, end_idx)
                chunk_similarities = similarity_values[start_idx:end_idx]
                
                rows = [f"{pair[0]},{pair[1]},{sim:.6f}" for pair, sim in zip(chunk_pairs, chunk_similarities)]
                f.write("\n".join(rows) + "\n")
                
                # Update the progress bar
                pbar.update(end_idx - start_idx)

"""Parallelize write out"""
def write_chunk(chunk_id, chunk_pairs, chunk_similarities, output_dir):
    """Write a single chunk to a temporary file."""
    temp_file = os.path.join(output_dir, f"temp_{chunk_id}.csv")
    with open(temp_file, 'w') as f:
        rows = [f"{pair[0]},{pair[1]},{sim:.6f}" for pair, sim in zip(chunk_pairs, chunk_similarities)]
        f.write("\n".join(rows) + "\n")
    return temp_file

def merge_temp_files(output_file, temp_files):
    """Merge all temporary files into the final output file."""
    with open(output_file, 'w') as final_file:
        # Write the header
        final_file.write("gene_1,gene_2,cosine_similarity\n")
        for temp_file in temp_files:
            with open(temp_file, 'r') as temp_f:
                final_file.write(temp_f.read())
            os.remove(temp_file)  # Clean up temporary file

def write_similarity_table_parallel(output_file, label_pairs, similarity_values, chunk_size=1000, num_workers=8):
    """Parallelize writing the similarity table."""
    assert len(label_pairs) == len(similarity_values), "Mismatch in label pairs and similarity values lengths."

    # create output directory for temporary files
    output_dir = os.path.dirname(output_file)
    if not output_dir:  # If no directory is specified, create tmp from the current directory
        output_dir = "./tmp/"
    os.makedirs(output_dir, exist_ok=True)

    # Split the data into chunks
    total_chunks = (len(label_pairs) + chunk_size - 1) // chunk_size
    chunks = [
        (chunk_id,
         label_pairs[chunk_id * chunk_size: (chunk_id + 1) * chunk_size],
         similarity_values[chunk_id * chunk_size: (chunk_id + 1) * chunk_size])
        for chunk_id in range(total_chunks)
    ]

    # Process chunks in parallel
    temp_files = []
    with Pool(processes=num_workers) as pool:
        with tqdm(total=len(chunks), desc="Writing chunks", unit="chunk") as pbar:
            for temp_file in pool.starmap(write_chunk, [(chunk_id, chunk_pairs, chunk_sims, output_dir)
                                                        for chunk_id, chunk_pairs, chunk_sims in chunks]):
                temp_files.append(temp_file)
                pbar.update(1)

    # Merge all temporary files into the final output
    merge_temp_files(output_file, temp_files)

def compute_all_pairs_cosine_similarity(embedding_dict, output_file, chunk_size=1000):
    start = time.time()
    labels = list(embedding_dict.keys())
    embeddings = torch.stack(list(embedding_dict.values()))
    embeddings = normalize(embeddings, p=2, dim=1)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    embeddings = embeddings.to(device)
    num_embeddings = embeddings.shape[0]
    similarity_matrix = []

    with tqdm(total=num_embeddings, desc="Computing similarity matrix", unit="rows") as pbar:
        for start_idx in range(0, num_embeddings, chunk_size):
            end_idx = min(start_idx + chunk_size, num_embeddings)
            chunk = embeddings[start_idx:end_idx]
            similarity_chunk = torch.mm(chunk, embeddings.T).cpu()  # Perform matrix multiplication on the GPU
            similarity_matrix.append(similarity_chunk)
            pbar.update(end_idx - start_idx)

    similarity_matrix = torch.cat(similarity_matrix, dim=0)

    print("Final similarity matrix created")
    print(similarity_matrix)

    print("Making label pairs")
    label_pairs = list(product(labels, repeat=2))
    print("Get similarity values")
    similarity_values = similarity_matrix.flatten().tolist()

    end = time.time()

    cos_sim_calc = end - start
    print('Cosine similarity calculation elapsed time:', cos_sim_calc)

    start_dict = time.time()

    # Create dictionary with label pairs as keys and similarity values as values
    similarity_dict = {
        (pair[0], pair[1]): sim for pair, sim in zip(label_pairs, similarity_values)
    }

    # Save dictionary as a .pkl file
    dict_output_file = output_file.replace('.csv', '_dict.pkl')  # Replace .csv with _dict.pkl
    with open(dict_output_file, 'wb') as f:
        pickle.dump(similarity_dict, f)

    print(f"Cosine similarity dictionary saved to {dict_output_file}")

    end_dict = time.time()

    make_dict = end_dict - start_dict
    print(f'Export to {dict_output_file} elapsed time:', make_dict)

    #print("Writing similarity table in parallel with progress tracking")
    # config for 208 core CPU available
    #write_similarity_table_parallel(output_file, label_pairs, similarity_values, chunk_size=100000, num_workers=208)
    #print(f"Cosine similarity table saved to {output_file}")

def create_synthetic_embeddings(output_pkl, num_embeddings=1000, embedding_dim=128):
    embedding_dict = {
        f"embedding_{i}": np.random.rand(embedding_dim) for i in range(num_embeddings)
    }
    with open(output_pkl, 'wb') as f:
        pickle.dump(embedding_dict, f)
    print(f"Synthetic embeddings saved to {output_pkl}")


def main():
    parser = argparse.ArgumentParser(description="Compute all-pairs cosine similarity from embedding vectors.")
    parser.add_argument('pkl_file', nargs='?', type=str, help="Path to the input .pkl file containing the embeddings dictionary.")
    parser.add_argument('output_file', nargs='?', type=str, help="Path to the output CSV file to save the results.")
    parser.add_argument('--create_synthetic', type=str, help="Create a synthetic embeddings file for testing.")
    parser.add_argument('--num_embeddings', type=int, default=1000, help="Number of synthetic embeddings (default: 1000).")
    parser.add_argument('--embedding_dim', type=int, default=128, help="Dimension of synthetic embeddings (default: 128).")
    parser.add_argument('--chunk_size', type=int, default=1000, help="Chunk size for matrix multiplication (default: 1000).")
    args = parser.parse_args()

    start_time = time.time()

    if args.create_synthetic:
        create_synthetic_embeddings(args.create_synthetic, args.num_embeddings, args.embedding_dim)
    elif args.pkl_file and args.output_file:
        embedding_dict = load_embeddings_from_pkl(args.pkl_file)
        compute_all_pairs_cosine_similarity(embedding_dict, args.output_file, args.chunk_size)
    else:
        parser.error("You must provide either --create_synthetic or both pkl_file and output_file.")

    end_time = time.time()  # Stop the timer
    elapsed_time = end_time - start_time
    print(f"Total script execution time: {elapsed_time:.2f} seconds")

if __name__ == '__main__':
    main()
