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
import h5py

def ensure_tensor(embedding_dict):
    """Ensure all embeddings in the dictionary are PyTorch tensors."""
    for key, value in embedding_dict.items():
        if isinstance(value, np.ndarray) or isinstance(value, list):
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

    print('Starting: write_similarity_table_parallel')
    
    assert len(label_pairs) == len(similarity_values), "Mismatch in label pairs and similarity values lengths."

    # # create output directory for temporary files
    # output_dir = os.path.dirname(output_file)
    # if not output_dir:  # If no directory is specified, create tmp from the current directory
    #     output_dir = "./tmp/"
    os.makedirs("./tmp/", exist_ok=True)

    print('output directory should be made')
    output_dir = "./tmp"

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

def compute_all_pairs_cosine_similarity(embedding_dict, output_file, chunk_size=1000, pairwise_df=None):
    start = time.time()
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    # Normalize embeddings and move to device
    labels = list(embedding_dict.keys())
    embeddings = torch.stack(list(embedding_dict.values()))
    embeddings = normalize(embeddings, p=2, dim=1)
    embeddings = embeddings.to(device)

    # Build lookup dictionary for embeddings
    embedding_lookup = {label: emb for label, emb in zip(labels, embeddings)}

    # Case 1: pairwise_df is provided
    if pairwise_df is not None:
        print("Computing cosine similarities for specific pairs from provided dataframe.")
        gene_a_list = pairwise_df['Gene_A'].tolist()
        gene_b_list = pairwise_df['Gene_B'].tolist()
    
        emb_a_list = []
        emb_b_list = []
        missing_genes = set()
    
        # Gather all embeddings into two lists
        for gene_a, gene_b in zip(gene_a_list, gene_b_list):
            emb_a = embedding_lookup.get(gene_a)
            emb_b = embedding_lookup.get(gene_b)
    
            if emb_a is None or emb_b is None:
                missing_genes.update([g for g in (gene_a, gene_b) if g not in embedding_lookup])
                emb_a_list.append(None)
                emb_b_list.append(None)
            else:
                emb_a_list.append(emb_a)
                emb_b_list.append(emb_b)
    
        # Indices for valid (non-missing) entries
        valid_indices = [i for i, (a, b) in enumerate(zip(emb_a_list, emb_b_list)) if a is not None and b is not None]
    
        # Prepare final list
        final_cosine_similarities = [float('nan')] * len(gene_a_list)
    
        batch_size = 100000  # or smaller if needed
    
        print(f"Calculating cosine similarities in batches of {batch_size}...")
    
        with tqdm(total=len(valid_indices), desc="Calculating batched cosine similarities", unit="pair") as pbar:
            for start in range(0, len(valid_indices), batch_size):
                end = min(start + batch_size, len(valid_indices))
                batch_indices = valid_indices[start:end]
    
                emb_a_batch = torch.stack([emb_a_list[i] for i in batch_indices]).to(device)
                emb_b_batch = torch.stack([emb_b_list[i] for i in batch_indices]).to(device)
    
                cos_batch = torch.sum(emb_a_batch * emb_b_batch, dim=1).cpu().tolist()
    
                # Write results into final list
                for idx, cos_sim in zip(batch_indices, cos_batch):
                    final_cosine_similarities[idx] = cos_sim
    
                pbar.update(len(batch_indices))
    
        if missing_genes:
            print(f"Warning: Missing embeddings for {len(missing_genes)} gene(s): {missing_genes}")
    
        # Build new dataframe
        result_df = pd.DataFrame({
            'Gene_A': gene_a_list,
            'Gene_B': gene_b_list,
            'Cosine_Similarity': final_cosine_similarities
        })
    
        # Save as pickle
        output_pkl = output_file.replace('.csv', '_cosine.pkl')
        result_df.to_pickle(output_pkl)
        print(f"Cosine similarity dataframe saved to {output_pkl}")

    # Case 2: default all-vs-all behavior
    else:
        print("Computing full all-vs-all cosine similarity matrix.")
        num_embeddings = embeddings.shape[0]
        similarity_matrix = []

        with tqdm(total=num_embeddings, desc="Computing similarity matrix", unit="rows") as pbar:
            for start_idx in range(0, num_embeddings, chunk_size):
                end_idx = min(start_idx + chunk_size, num_embeddings)
                chunk = embeddings[start_idx:end_idx]
                similarity_chunk = torch.mm(chunk, embeddings.T).cpu()
                similarity_matrix.append(similarity_chunk)
                pbar.update(end_idx - start_idx)

        similarity_matrix = torch.cat(similarity_matrix, dim=0)

        print("Final similarity matrix created.")

        label_pairs = list(product(labels, repeat=2))
        similarity_values = similarity_matrix.flatten().tolist()

        dict_output_file = output_file.replace('.csv', '_dict.h5')
        chunk_size_h5 = min(1000000, len(label_pairs))

        with h5py.File(dict_output_file, 'w') as f:
            dset_pairs = f.create_dataset("gene_pairs", (len(label_pairs), 2), dtype="S50", chunks=True)
            dset_similarities = f.create_dataset("similarities", (len(label_pairs),), dtype="float32", chunks=True)

            for i in range(0, len(label_pairs), chunk_size_h5):
                chunk_pairs = np.array(label_pairs[i:i + chunk_size_h5], dtype="S50")
                chunk_similarities = np.array(similarity_values[i:i + chunk_size_h5], dtype="float32")
                dset_pairs[i:i + chunk_size_h5] = chunk_pairs
                dset_similarities[i:i + chunk_size_h5] = chunk_similarities
                print(f"Processed chunk {i // chunk_size_h5 + 1}/{(len(label_pairs) + chunk_size_h5 - 1) // chunk_size_h5}")

        print(f"Cosine similarity dictionary saved to {dict_output_file}")

    end = time.time()
    print(f"Total cosine similarity calculation elapsed time: {end - start:.2f} seconds")

def create_synthetic_embeddings(output_pkl, num_embeddings=1000, embedding_dim=128):
    embedding_dict = {
        f"embedding_{i}": np.random.rand(embedding_dim) for i in range(num_embeddings)
    }
    with open(output_pkl, 'wb') as f:
        pickle.dump(embedding_dict, f)
    print(f"Synthetic embeddings saved to {output_pkl}")

def main():
    parser = argparse.ArgumentParser(description="Compute all-pairs cosine similarity from embedding vectors.")
    parser.add_argument('embedding_pkl', type=str, help="Path to the input .pkl file containing the embeddings dictionary.")
    parser.add_argument('output_file', type=str, help="Path to the output file to save the results (format determined automatically).")
    parser.add_argument('--pairwise_pkl', type=str, default=None, help="Optional: Path to a .pkl file containing specific gene pairs to compute.")
    parser.add_argument('--chunk_size', type=int, default=1000, help="Chunk size for matrix multiplication (default: 1000).")
    parser.add_argument('--create_synthetic', type=str, help="Optional: Path to create a synthetic embeddings file for testing.")
    parser.add_argument('--num_embeddings', type=int, default=1000, help="Number of synthetic embeddings (default: 1000).")
    parser.add_argument('--embedding_dim', type=int, default=128, help="Dimension of synthetic embeddings (default: 128).")
    args = parser.parse_args()

    start_time = time.time()

    if args.create_synthetic:
        create_synthetic_embeddings(args.create_synthetic, args.num_embeddings, args.embedding_dim)
    else:
        embedding_dict = load_embeddings_from_pkl(args.embedding_pkl)

        pairwise_df = None
        if args.pairwise_pkl:
            print(f"Loading pairwise gene list from {args.pairwise_pkl}")
            with open(args.pairwise_pkl, 'rb') as f:
                pairwise_df = pickle.load(f)
            if not isinstance(pairwise_df, pd.DataFrame):
                raise ValueError("The pairwise_pkl file must contain a pandas DataFrame with Gene_A and Gene_B columns.")

        compute_all_pairs_cosine_similarity(embedding_dict, args.output_file, args.chunk_size, pairwise_df=pairwise_df)

    end_time = time.time()
    print(f"Total script execution time: {end_time - start_time:.2f} seconds")

if __name__ == '__main__':
    main()
